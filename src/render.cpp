#include "raym0nade/render.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "raym0nade/image.hpp"
#include "raym0nade/sampling.hpp"

namespace raym0nade {
namespace {

constexpr int kMaximumPathDepth = 16;
constexpr int kMaximumSamplesPerPixel = 1'000'000;
constexpr int kTransparentSampleMultiplier = 16;
constexpr float kRegularizationFactor = 1.0F;
constexpr float kBaseReflectionProbability = 0.24F;
constexpr float kMaximumExposure = 1.0e6F;

struct MediumEntry {
    int materialId{-1};
    float ior{1.0F};
    vec3 absorption{1.0F};
};

class MediumStack {
public:
    void reset() {
        entries_.clear();
        entries_.push_back(MediumEntry{});
    }

    void enter(int materialId, float ior, const vec3& absorption) {
        entries_.push_back(MediumEntry{materialId, std::max(ior, 1.0e-4F), absorption});
    }

    void exit(int materialId) {
        for (auto entry = entries_.rbegin(); entry != entries_.rend(); ++entry) {
            if (entry->materialId == materialId) {
                entries_.erase(std::next(entry).base());
                return;
            }
        }
    }

    [[nodiscard]] float currentIor() const noexcept {
        return entries_.empty() ? 1.0F : entries_.back().ior;
    }

    [[nodiscard]] float iorAfterExit(int materialId) const noexcept {
        for (auto entry = entries_.rbegin(); entry != entries_.rend(); ++entry) {
            if (entry->materialId == materialId) {
                const auto forward = std::next(entry).base();
                return forward == entries_.begin() ? 1.0F : std::prev(forward)->ior;
            }
        }
        return 1.0F;
    }

    [[nodiscard]] vec3 combinedAbsorption() const noexcept {
        vec3 value{1.0F};
        for (const MediumEntry& entry : entries_) {
            value *= entry.absorption;
        }
        return value;
    }

    [[nodiscard]] std::size_t size() const noexcept {
        return entries_.size();
    }

private:
    std::vector<MediumEntry> entries_{MediumEntry{}};
};

struct RenderContext {
    RenderContext(
        std::uint32_t seed,
        DirectLightSamplingScratch& primaryDirectLightSamplingScratch,
        DirectLightSamplingScratch& pathDirectLightSamplingScratch)
        : generator(seed),
          primaryDirectLightSamplingScratch(primaryDirectLightSamplingScratch),
          pathDirectLightSamplingScratch(pathDirectLightSamplingScratch) {
        lightSamples.reserve(6);
    }

    Generator generator;
    MediumStack media;
    std::vector<LightSample> lightSamples;
    DirectLightSamplingScratch& primaryDirectLightSamplingScratch;
    DirectLightSamplingScratch& pathDirectLightSamplingScratch;
    std::uint64_t directLightSamples{0};
};

std::uint32_t hashSeed(std::uint32_t seed, std::uint32_t row) noexcept {
    std::uint32_t value = seed ^ (row + 0x9E3779B9U + (seed << 6U) + (seed >> 2U));
    value ^= value >> 16U;
    value *= 0x7FEB352DU;
    value ^= value >> 15U;
    value *= 0x846CA68BU;
    value ^= value >> 16U;
    return value;
}

void calculatePositionDifferentials(
    const Ray& ray,
    float hitDistance,
    const vec3& normal,
    const RayDifferential& base,
    vec3& positionDx,
    vec3& positionDy) noexcept {
    const float denominator = glm::dot(ray.direction, normal);
    if (!isFinite(denominator) || std::abs(denominator) <= 1.0e-8F) {
        positionDx = vec3{0.0F};
        positionDy = vec3{0.0F};
        return;
    }
    const float distanceDx =
        -glm::dot(base.dPdx + hitDistance * base.dDdx, normal) / denominator;
    const float distanceDy =
        -glm::dot(base.dPdy + hitDistance * base.dDdy, normal) / denominator;
    positionDx = base.dPdx + distanceDx * ray.direction + hitDistance * base.dDdx;
    positionDy = base.dPdy + distanceDy * ray.direction + hitDistance * base.dDdy;
}

void calculateDirectionDifferentials(
    const Ray& ray,
    const vec3& normal,
    const RayDifferential& base,
    vec3& directionDx,
    vec3& directionDy) noexcept {
    const float normalDx = glm::dot(base.dDdx, normal);
    const float normalDy = glm::dot(base.dDdy, normal);
    directionDx = base.dDdx - 2.0F * normalDx * normal;
    directionDy = base.dDdy - 2.0F * normalDy * normal;
    if (!isFinite(directionDx)) {
        directionDx = vec3{0.0F};
    }
    if (!isFinite(directionDy)) {
        directionDy = vec3{0.0F};
    }
    (void)ray;
}

void populateHitInfo(
    const HitRecord& hit,
    const Ray& ray,
    const RayDifferential& differential,
    vec3& positionDx,
    vec3& positionDy,
    HitInfo& hitInfo) {
    const Face& face = *hit.face;
    const vec3 coordinates =
        barycentric(face.vertices[0], face.vertices[1], face.vertices[2], hitInfo.position);
    getHitNormals(
        face, ray.direction, coordinates, hitInfo.shapeNormal, hitInfo.surfaceNormal, hitInfo.entering);
    const vec3 originalSurfaceNormal = hitInfo.surfaceNormal;
    calculatePositionDifferentials(
        ray, hit.tMaximum, hitInfo.shapeNormal, differential, positionDx, positionDy);
    getHitMaterial(face, coordinates, positionDx, positionDy, hitInfo);
    if (glm::dot(hitInfo.surfaceNormal, ray.direction) >= 0.0F) {
        hitInfo.surfaceNormal = originalSurfaceNormal;
    }
    if (glm::dot(hitInfo.surfaceNormal, ray.direction) >= 0.0F) {
        hitInfo.surfaceNormal = hitInfo.shapeNormal;
    }
}

vec3 mediumTransmission(const vec3& absorption, float distance) noexcept {
    constexpr float absorptionScale = 32.0F;
    if (!isFinite(absorption) || !isFinite(distance)) {
        return vec3{1.0F};
    }
    const vec3 boundedAbsorption = glm::clamp(absorption, vec3{0.0F}, vec3{1.0F});
    if (boundedAbsorption.x >= 1.0F - kRayEpsilon &&
        boundedAbsorption.y >= 1.0F - kRayEpsilon &&
        boundedAbsorption.z >= 1.0F - kRayEpsilon) {
        return vec3{1.0F};
    }
    const double exponent = static_cast<double>(absorptionScale) *
                            static_cast<double>(std::max(distance, 0.0F));
    const float boundedExponent = static_cast<float>(std::min(
        exponent, static_cast<double>(std::numeric_limits<float>::max())));
    return positivePow(boundedAbsorption, boundedExponent);
}

float relativeEta(const MediumStack& media, const HitInfo& hit) noexcept {
    const float from = media.currentIor();
    const float to = hit.entering ? hit.eta : media.iorAfterExit(hit.materialId);
    return from / std::max(to, 1.0e-4F);
}

void applyTransmissionBoundary(MediumStack& media, const HitInfo& hit, float absoluteIor) {
    if (hit.entering) {
        media.enter(hit.materialId, absoluteIor, hit.baseColor);
    } else {
        media.exit(hit.materialId);
    }
}

void scaleThroughput(std::vector<LightSample>& samples, float factor) noexcept {
    for (LightSample& sample : samples) {
        sample.throughput *= factor;
    }
}

void propagateThroughput(std::vector<LightSample>& samples, const vec3& throughput) noexcept {
    for (LightSample& sample : samples) {
        sample.radiance *= sample.throughput;
        sample.throughput = throughput;
    }
}

void scaleWeights(std::vector<LightSample>& samples, float factor) noexcept {
    for (LightSample& sample : samples) {
        sample.weight *= factor;
    }
}

void traceRay(
    const Ray& ray,
    const RayDifferential& differential,
    const Model& model,
    RenderContext& context,
    float roughnessFloor,
    bool excludeDirectLight,
    int depth,
    std::vector<LightSample>& samples) {
    static constexpr std::array<int, kMaximumPathDepth + 1> lightSampleCounts{
        0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6};

    const HitRecord hit = model.intersect(ray);
    if (hit.face == nullptr) {
        if (excludeDirectLight || model.sky().empty()) {
            return;
        }
        samples.push_back(
            LightSample{vec3{1.0F}, model.sky().radiance(ray.direction), 1.0F});
        return;
    }

    Bsdf bsdf{-ray.direction, ray.origin + ray.direction * hit.tMaximum};
    vec3 positionDx;
    vec3 positionDy;
    populateHitInfo(hit, ray, differential, positionDx, positionDy, bsdf.surface);
    const float absoluteIor = bsdf.surface.eta;

    roughnessFloor = std::max(roughnessFloor, kRegularizationFactor * bsdf.surface.roughness);
    bsdf.surface.roughness = std::max(bsdf.surface.roughness, roughnessFloor);
    const float rouletteProbability = std::clamp(
        1.0F - 0.5F * safeSqrt(bsdf.surface.roughness), 0.5F, 1.0F);
    bool maySampleDirect = rouletteProbability < 0.9F;
    const vec3 absorption = mediumTransmission(context.media.combinedAbsorption(), hit.tMaximum);

    if (glm::length(bsdf.surface.emission) > kRayEpsilon) {
        if (excludeDirectLight) {
            return;
        }
        samples.push_back(LightSample{absorption, bsdf.surface.emission, 1.0F});
        return;
    }

    float reflectionProbability = 1.0F;
    float fresnel = 1.0F;
    vec3 exactRefraction{std::numeric_limits<float>::quiet_NaN()};
    if (bsdf.surface.opacity < kRayEpsilon) {
        bsdf.surface.eta = relativeEta(context.media, bsdf.surface);
        bsdf.refractPrecisely(exactRefraction, fresnel);
        reflectionProbability = context.media.size() == 1
                                    ? kBaseReflectionProbability +
                                          (1.0F - kBaseReflectionProbability) * fresnel
                                    : fresnel;
        reflectionProbability = std::clamp(reflectionProbability, 0.0F, 1.0F);
    }

    vec3 newDirection{std::numeric_limits<float>::quiet_NaN()};
    vec3 throughput{0.0F};
    RayDifferential nextDifferential;
    int failures = 0;
    const bool reflect = reflectionProbability >= 1.0F ||
                         (reflectionProbability > 0.0F && context.generator() < reflectionProbability);

    if (reflect) {
        maySampleDirect = maySampleDirect && bsdf.surface.entering && context.media.size() == 1;
        if (maySampleDirect &&
            (depth >= kMaximumPathDepth || context.generator() > rouletteProbability)) {
            ++context.directLightSamples;
            sampleDirectLight(
                bsdf,
                model,
                context.generator,
                lightSampleCounts[static_cast<std::size_t>(depth)],
                samples,
                context.pathDirectLightSamplingScratch);
            float scale = 1.0F / std::max(reflectionProbability, 1.0e-6F);
            if (depth < kMaximumPathDepth) {
                scale /= std::max(1.0F - rouletteProbability, 1.0e-6F);
            }
            scaleThroughput(samples, scale);
            return;
        }

        bsdf.sampleReflection(context.generator, newDirection, throughput, failures);
        throughput /= std::max(reflectionProbability, 1.0e-6F);
        if (maySampleDirect) {
            throughput /= std::max(rouletteProbability, 1.0e-6F);
        }
        vec3 directionDx;
        vec3 directionDy;
        calculateDirectionDifferentials(
            ray, bsdf.surface.surfaceNormal, differential, directionDx, directionDy);
        nextDifferential = RayDifferential{positionDx, positionDy, directionDx, directionDy};
    } else {
        const float transmissionProbability = 1.0F - reflectionProbability;
        maySampleDirect = maySampleDirect && !bsdf.surface.entering && context.media.size() == 2;
        if (maySampleDirect &&
            (depth >= kMaximumPathDepth || context.generator() > rouletteProbability)) {
            ++context.directLightSamples;
            sampleDirectLight(
                bsdf,
                model,
                context.generator,
                lightSampleCounts[static_cast<std::size_t>(depth)],
                samples,
                context.pathDirectLightSamplingScratch);
            float scale = 1.0F / std::max(transmissionProbability, 1.0e-6F);
            if (depth < kMaximumPathDepth) {
                scale /= std::max(1.0F - rouletteProbability, 1.0e-6F);
            }
            scaleThroughput(samples, scale);
            propagateThroughput(samples, absorption);
            return;
        }

        if (std::abs(bsdf.surface.eta - 1.0F) < kRayEpsilon) {
            newDirection = ray.direction;
            throughput = vec3{1.0F};
            nextDifferential = differential;
        } else {
            bsdf.sampleTransmission(context.generator, newDirection, throughput, failures);
            nextDifferential = differential;
        }
        throughput /= std::max(transmissionProbability, 1.0e-6F);
        if (maySampleDirect) {
            throughput /= std::max(rouletteProbability, 1.0e-6F);
        }
        applyTransmissionBoundary(context.media, bsdf.surface, absoluteIor);
    }

    (void)failures;
    if (!isFinite(newDirection) || !isFinite(throughput) || depth >= kMaximumPathDepth) {
        return;
    }
    throughput *= absorption;
    traceRay(
        Ray{bsdf.surface.position, newDirection},
        nextDifferential,
        model,
        context,
        roughnessFloor,
        maySampleDirect,
        depth + 1,
        samples);
    propagateThroughput(samples, throughput);
}

void sampleIndirectFromFirstHit(
    const HitInfo& hitInfo,
    const vec3& origin,
    const RayDifferential& differential,
    const Model& model,
    RenderContext& context,
    std::vector<LightSample>& samples) {
    samples.clear();
    if (glm::length(hitInfo.emission) > 0.0F) {
        return;
    }
    context.media.reset();

    const vec3 cameraOffset = hitInfo.position - origin;
    const float hitDistance = glm::length(cameraOffset);
    const vec3 incomingDirection = safeNormalize(cameraOffset);
    if (glm::dot(incomingDirection, incomingDirection) == 0.0F) {
        return;
    }
    const Ray incomingRay{origin, incomingDirection};
    Bsdf bsdf{-incomingDirection, hitInfo};
    const float roughnessFloor = hitInfo.roughness * kRegularizationFactor;
    const float absoluteIor = bsdf.surface.eta;
    bsdf.surface.eta = 1.0F / std::max(bsdf.surface.eta, 1.0e-4F);

    float reflectionProbability = 1.0F;
    float fresnel = 1.0F;
    vec3 exactRefraction;
    if (bsdf.surface.opacity < kRayEpsilon) {
        bsdf.refractPrecisely(exactRefraction, fresnel);
        reflectionProbability = std::clamp(
            kBaseReflectionProbability + (1.0F - kBaseReflectionProbability) * fresnel,
            0.0F,
            1.0F);
    }

    vec3 newDirection{std::numeric_limits<float>::quiet_NaN()};
    vec3 throughput{0.0F};
    RayDifferential nextDifferential;
    int failures = 0;
    const bool reflect = reflectionProbability >= 1.0F ||
                         (reflectionProbability > 0.0F && context.generator() < reflectionProbability);
    if (reflect) {
        vec3 positionDx;
        vec3 positionDy;
        calculatePositionDifferentials(
            incomingRay, hitDistance, bsdf.surface.shapeNormal, differential, positionDx, positionDy);
        vec3 directionDx;
        vec3 directionDy;
        calculateDirectionDifferentials(
            incomingRay, bsdf.surface.surfaceNormal, differential, directionDx, directionDy);
        nextDifferential = RayDifferential{positionDx, positionDy, directionDx, directionDy};
        bsdf.sampleReflection(context.generator, newDirection, throughput, failures);
        throughput /= std::max(reflectionProbability, 1.0e-6F);
    } else {
        bsdf.sampleTransmission(context.generator, newDirection, throughput, failures);
        throughput /= std::max(1.0F - reflectionProbability, 1.0e-6F);
        nextDifferential = differential;
        applyTransmissionBoundary(context.media, bsdf.surface, absoluteIor);
    }

    (void)failures;
    if (!isFinite(newDirection) || !isFinite(throughput)) {
        return;
    }
    traceRay(
        Ray{bsdf.surface.position, newDirection},
        nextDifferential,
        model,
        context,
        roughnessFloor,
        true,
        1,
        samples);
    propagateThroughput(samples, throughput);
}

void sampleDirectFromFirstHit(
    const HitInfo& hitInfo,
    const vec3& origin,
    const Model& model,
    int sampleCount,
    float techniqueScale,
    float sampleWeight,
    RenderContext& context,
    RadianceData& diffuseRadiance,
    RadianceData& specularRadiance) {
    if (!isFinite(hitInfo.position) || glm::length(hitInfo.emission) > 0.0F || sampleCount <= 0 ||
        !isFinite(techniqueScale) || techniqueScale <= 0.0F || !isFinite(sampleWeight) ||
        sampleWeight <= 0.0F) {
        return;
    }
    const vec3 incoming = safeNormalize(hitInfo.position - origin);
    if (glm::dot(incoming, incoming) == 0.0F) {
        return;
    }
    HitInfo directHit = hitInfo;
    if (directHit.opacity < kRayEpsilon) {
        directHit.eta = 1.0F / std::max(directHit.eta, 1.0e-4F);
    }
    const Bsdf bsdf{-incoming, directHit};
    sampleDirectLight(
        bsdf,
        model,
        context.generator,
        sampleCount,
        context.lightSamples,
        context.primaryDirectLightSamplingScratch);
    for (LightSample sample : context.lightSamples) {
        sample.throughput *= techniqueScale;
        sample.weight *= sampleWeight;
        accumulateInwardRadiance(hitInfo.baseColor, sample, diffuseRadiance, specularRadiance);
    }
}

RayDifferential initialRayDifferential(
    const vec3& unnormalizedDirection, const PinholeCamera& camera) noexcept {
    RayDifferential differential;
    const vec3 directionDx = camera.pixelScale * camera.right;
    const vec3 directionDy = camera.pixelScale * camera.up;
    const float squaredLength = glm::dot(unnormalizedDirection, unnormalizedDirection);
    if (squaredLength <= 1.0e-12F || !isFinite(squaredLength)) {
        return differential;
    }
    const float denominator = safeSqrt(squaredLength) * squaredLength;
    differential.dDdx =
        (squaredLength * directionDx - unnormalizedDirection * glm::dot(unnormalizedDirection, directionDx)) /
        denominator;
    differential.dDdy =
        (squaredLength * directionDy - unnormalizedDirection * glm::dot(unnormalizedDirection, directionDy)) /
        denominator;
    return differential;
}

void renderPixel(
    const Model& model,
    const RenderSettings& settings,
    const PinholeCamera& camera,
    const ImageExtent& extent,
    RenderContext& context,
    Film& film,
    int x,
    int y) {
    const std::size_t pixelIndex = static_cast<std::size_t>(y) * static_cast<std::size_t>(settings.width) +
                                   static_cast<std::size_t>(x);
    const auto pixelX = static_cast<std::uint32_t>(x);
    const auto pixelY = static_cast<std::uint32_t>(y);
    const vec3 unnormalizedDirection =
        primaryRayDirectionUnnormalized(camera, extent, pixelX, pixelY);
    const RayDifferential differential = initialRayDifferential(unnormalizedDirection, camera);
    const Ray ray{camera.position, safeNormalize(unnormalizedDirection)};
    HitInfo& gBuffer = film.gBuffer[pixelIndex];

    const HitRecord hit = model.intersect(ray);
    if (hit.face == nullptr) {
        gBuffer.position = vec3{std::numeric_limits<float>::quiet_NaN()};
        if (!model.sky().empty()) {
            gBuffer.emission = model.sky().radiance(ray.direction);
        }
        return;
    }

    gBuffer.position = ray.origin + hit.tMaximum * ray.direction;
    vec3 positionDx;
    vec3 positionDy;
    populateHitInfo(hit, ray, differential, positionDx, positionDy, gBuffer);

    const bool opaque = gBuffer.opacity > 1.0F - kRayEpsilon;
    RadianceData& directDiffuse = film.directDiffuseRadiance[pixelIndex];
    RadianceData& directSpecular = film.directSpecularRadiance[pixelIndex];
    RadianceData& indirectDiffuse = film.indirectDiffuseRadiance[pixelIndex];
    RadianceData& indirectSpecular = film.indirectSpecularRadiance[pixelIndex];

    const float directProbability = settings.directLightProbability;
    const float indirectProbability = 1.0F - directProbability;
    const float sampleWeight = 1.0F / static_cast<float>(settings.samplesPerPixel);
    const int indirectReplicates = opaque ? 1 : kTransparentSampleMultiplier;
    for (int sampleIndex = 0; sampleIndex < settings.samplesPerPixel; ++sampleIndex) {
        const bool chooseDirect = directProbability >= 1.0F ||
                                  (directProbability > 0.0F &&
                                   context.generator() < directProbability);
        if (chooseDirect) {
            ++context.directLightSamples;
            sampleDirectFromFirstHit(
                gBuffer,
                settings.position,
                model,
                1,
                1.0F / directProbability,
                sampleWeight,
                context,
                directDiffuse,
                directSpecular);
            continue;
        }

        for (int replicate = 0; replicate < indirectReplicates; ++replicate) {
            sampleIndirectFromFirstHit(
                gBuffer,
                settings.position,
                differential,
                model,
                context,
                context.lightSamples);
            scaleThroughput(context.lightSamples, 1.0F / indirectProbability);
            scaleWeights(
                context.lightSamples, sampleWeight / static_cast<float>(indirectReplicates));
            for (const LightSample& sample : context.lightSamples) {
                accumulateInwardRadiance(
                    gBuffer.opacity < kRayEpsilon ? vec3{0.0F} : gBuffer.baseColor,
                    sample,
                    indirectDiffuse,
                    indirectSpecular);
            }
        }
    }

    finalizeRadianceData(directDiffuse, settings.exposure);
    finalizeRadianceData(directSpecular, settings.exposure);
    finalizeRadianceData(indirectDiffuse, settings.exposure);
    finalizeRadianceData(indirectSpecular, settings.exposure);
}

void renderRows(
    const Model& model,
    const RenderSettings& settings,
    Film& film,
    std::atomic<int>& nextRow,
    std::atomic<std::uint64_t>& directLightSamples) {
    const PinholeCamera camera{
        settings.position,
        settings.direction,
        settings.up,
        settings.right,
        settings.pixelScale,
    };
    const ImageExtent extent{
        static_cast<std::uint32_t>(settings.width),
        static_cast<std::uint32_t>(settings.height),
    };
    DirectLightSamplingScratch primaryDirectLightSamplingScratch;
    DirectLightSamplingScratch pathDirectLightSamplingScratch;
    while (true) {
        const int row = nextRow.fetch_add(1, std::memory_order_relaxed);
        if (row >= settings.height) {
            return;
        }
        RenderContext context{
            hashSeed(settings.seed, static_cast<std::uint32_t>(row)),
            primaryDirectLightSamplingScratch,
            pathDirectLightSamplingScratch};
        for (int x = 0; x < settings.width; ++x) {
            renderPixel(model, settings, camera, extent, context, film, x, row);
        }
        directLightSamples.fetch_add(context.directLightSamples, std::memory_order_relaxed);
    }
}

FilmRenderResult renderToFilmValidated(
    const Model& model, const RenderSettings& settings, int workerCount) {
    Film film{settings.width, settings.height};
    film.exposure = settings.exposure;

    std::atomic<int> nextRow{0};
    std::atomic<std::uint64_t> directLightSamples{0};
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(workerCount));
    std::vector<std::exception_ptr> workerErrors(static_cast<std::size_t>(workerCount));
    const auto renderStart = std::chrono::steady_clock::now();
    try {
        for (int index = 0; index < workerCount; ++index) {
            workers.emplace_back([&, index] {
                try {
                    renderRows(model, settings, film, nextRow, directLightSamples);
                } catch (...) {
                    workerErrors[static_cast<std::size_t>(index)] = std::current_exception();
                }
            });
        }
    } catch (...) {
        nextRow.store(settings.height, std::memory_order_relaxed);
        for (std::thread& worker : workers) {
            worker.join();
        }
        throw;
    }
    for (std::thread& worker : workers) {
        worker.join();
    }
    for (const std::exception_ptr& error : workerErrors) {
        if (error != nullptr) {
            std::rethrow_exception(error);
        }
    }
    const auto renderEnd = std::chrono::steady_clock::now();

    RenderStats stats;
    stats.renderSeconds = std::chrono::duration<double>(renderEnd - renderStart).count();
    stats.totalSeconds = stats.renderSeconds;
    stats.directLightSamples = directLightSamples.load(std::memory_order_relaxed);
    return FilmRenderResult{std::move(film), stats};
}

void renderPrimaryAovRow(
    const Model& model,
    const PrimaryRenderRequest& request,
    const vec3& directionToLight,
    LinearImage& image,
    std::uint32_t y) {
    for (std::uint32_t x = 0U; x < request.extent.width; ++x) {
        const std::size_t pixelIndex =
            static_cast<std::size_t>(y) * static_cast<std::size_t>(request.extent.width) +
            static_cast<std::size_t>(x);
        const vec3 unnormalizedDirection =
            primaryRayDirectionUnnormalized(request.camera, request.extent, x, y);
        const Ray ray{request.camera.position, safeNormalize(unnormalizedDirection)};
        const HitRecord hit = model.intersect(ray);
        if (hit.face == nullptr) {
            continue;
        }

        HitInfo hitInfo;
        hitInfo.position = ray.origin + hit.tMaximum * ray.direction;
        if (request.aov == PrimaryAov::ShapeNormal) {
            const Face& face = *hit.face;
            const vec3 coordinates = barycentric(
                face.vertices[0], face.vertices[1], face.vertices[2], hitInfo.position);
            getHitNormals(
                face,
                ray.direction,
                coordinates,
                hitInfo.shapeNormal,
                hitInfo.surfaceNormal,
                hitInfo.entering);
            image.pixels[pixelIndex] = isFinite(hitInfo.shapeNormal)
                                           ? (hitInfo.shapeNormal + vec3{1.0F}) * 0.5F
                                           : vec3{0.0F};
            continue;
        }

        const RayDifferential differential =
            initialRayDifferential(unnormalizedDirection, request.camera);
        vec3 positionDx;
        vec3 positionDy;
        populateHitInfo(hit, ray, differential, positionDx, positionDy, hitInfo);

        if (request.aov == PrimaryAov::BaseColor) {
            image.pixels[pixelIndex] =
                isFinite(hitInfo.baseColor) ? hitInfo.baseColor : vec3{0.0F};
        } else {
            if (!isFinite(hitInfo.shapeNormal) || !isFinite(hitInfo.baseColor)) {
                continue;
            }
            const float normalDotLight = glm::dot(hitInfo.shapeNormal, directionToLight);
            if (!(normalDotLight > 0.0F)) {
                continue;
            }
            if (model.occluded(
                    Ray{hitInfo.position, directionToLight},
                    std::numeric_limits<float>::infinity())) {
                continue;
            }
            const vec3 radiance =
                glm::max(hitInfo.baseColor, vec3{0.0F}) *
                request.directionalLight.incidentRadiance * (normalDotLight / kPi);
            image.pixels[pixelIndex] = isFinite(radiance) ? radiance : vec3{0.0F};
        }
    }
}

void renderPrimaryAovRows(
    const Model& model,
    const PrimaryRenderRequest& request,
    const vec3& directionToLight,
    LinearImage& image,
    std::atomic<std::uint32_t>& nextRow) {
    while (true) {
        const std::uint32_t row = nextRow.fetch_add(1U, std::memory_order_relaxed);
        if (row >= request.extent.height) {
            return;
        }
        renderPrimaryAovRow(model, request, directionToLight, image, row);
    }
}

}  // namespace

LinearImage renderPrimaryAovCpu(const Model& model, const PrimaryRenderRequest& request) {
    return renderPrimaryAovCpu(model, request, CpuPrimaryRenderOptions{});
}

void CpuPrimaryRenderOptions::validate() const {
    if (threadCount < 0) {
        throw std::invalid_argument(
            "CPU primary-render thread count must be zero (automatic) or positive.");
    }
}

std::uint32_t CpuPrimaryRenderOptions::resolvedThreadCount(
    std::uint32_t imageHeight) const noexcept {
    const unsigned int hardwareThreads =
        std::max(1U, std::thread::hardware_concurrency());
    const std::uint32_t available = static_cast<std::uint32_t>(hardwareThreads);
    const std::uint32_t requested =
        threadCount > 0 ? static_cast<std::uint32_t>(threadCount) : available;
    const std::uint32_t boundedHeight = imageHeight == 0U ? 1U : imageHeight;
    return std::min(requested, boundedHeight);
}

LinearImage renderPrimaryAovCpu(
    const Model& model,
    const PrimaryRenderRequest& request,
    const CpuPrimaryRenderOptions& options) {
    request.validate();
    options.validate();
    const vec3 directionToLight = safeNormalize(request.directionalLight.directionToLight);
    LinearImage image{
        request.extent,
        std::vector<vec3>(request.extent.pixelCount(), vec3{0.0F}),
    };

    const std::uint32_t workerCount =
        options.resolvedThreadCount(request.extent.height);
    std::atomic<std::uint32_t> nextRow{0U};
    std::vector<std::thread> workers;
    workers.reserve(workerCount);
    std::vector<std::exception_ptr> workerErrors(workerCount);
    try {
        for (std::uint32_t index = 0U; index < workerCount; ++index) {
            workers.emplace_back([&, index] {
                try {
                    renderPrimaryAovRows(
                        model, request, directionToLight, image, nextRow);
                } catch (...) {
                    nextRow.store(request.extent.height, std::memory_order_relaxed);
                    workerErrors[index] = std::current_exception();
                }
            });
        }
    } catch (...) {
        nextRow.store(request.extent.height, std::memory_order_relaxed);
        for (std::thread& worker : workers) {
            worker.join();
        }
        throw;
    }
    for (std::thread& worker : workers) {
        worker.join();
    }
    for (const std::exception_ptr& error : workerErrors) {
        if (error != nullptr) {
            std::rethrow_exception(error);
        }
    }
    return image;
}

void RenderSettings::validate() const {
    if (width <= 0 || height <= 0) {
        throw std::invalid_argument("Render width and height must be positive.");
    }
    constexpr std::uint64_t maximumPixels = 1ULL << 30U;
    if (static_cast<std::uint64_t>(width) * static_cast<std::uint64_t>(height) > maximumPixels) {
        throw std::invalid_argument("Render dimensions are too large.");
    }
    if (samplesPerPixel <= 0) {
        throw std::invalid_argument("Samples per pixel must be positive.");
    }
    if (samplesPerPixel > kMaximumSamplesPerPixel) {
        throw std::invalid_argument("Samples per pixel exceed the supported limit.");
    }
    if (threadCount < 0) {
        throw std::invalid_argument("Thread count must be zero (automatic) or positive.");
    }
    if (!isFinite(position) || !isFinite(direction) || !isFinite(up) || !isFinite(right) ||
        glm::dot(direction, direction) <= 1.0e-12F || glm::dot(up, up) <= 1.0e-12F ||
        glm::dot(right, right) <= 1.0e-12F) {
        throw std::invalid_argument("Camera vectors must be finite and non-zero.");
    }
    const vec3 normalizedDirection = safeNormalize(direction);
    const vec3 normalizedRight = safeNormalize(right);
    const vec3 normalizedUp = safeNormalize(up);
    const float basisVolume = std::abs(glm::dot(
        glm::cross(normalizedDirection, normalizedRight), normalizedUp));
    if (!isFinite(basisVolume) || basisVolume <= 1.0e-4F) {
        throw std::invalid_argument("Camera direction, right, and up vectors must be linearly independent.");
    }
    if (!isFinite(pixelScale) || pixelScale <= 0.0F || !isFinite(exposure) || exposure < 0.0F ||
        exposure > kMaximumExposure || !isFinite(focusDistance) || focusDistance < 0.0F ||
        !isFinite(circleOfConfusion) || circleOfConfusion < 0.0F) {
        throw std::invalid_argument("Camera and exposure scalars are invalid.");
    }
    if (!isFinite(directLightProbability) || directLightProbability < 0.0F ||
        directLightProbability > 1.0F) {
        throw std::invalid_argument("Direct-light probability must be in the [0, 1] range.");
    }
    if (outputPrefix.empty()) {
        throw std::invalid_argument("Output prefix must not be empty.");
    }
}

int RenderSettings::resolvedThreadCount() const noexcept {
    const unsigned int available = std::max(1U, std::thread::hardware_concurrency());
    const int requested = threadCount == 0 ? static_cast<int>(available) : threadCount;
    return std::max(1, std::min(requested, height));
}

FilmRenderResult renderToFilm(const Model& model, const RenderSettings& settings) {
    settings.validate();
    return renderToFilmValidated(model, settings, settings.resolvedThreadCount());
}

void exportFilmToFiles(Film film, const RenderSettings& settings) {
    settings.validate();
    if (film.width() != settings.width || film.height() != settings.height) {
        throw std::invalid_argument(
            "Film dimensions must match the render settings before export.");
    }
    const std::filesystem::path parent = settings.outputPrefix.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }

    const auto saveCurrentImage = [&](const std::string& tag) {
        std::filesystem::path filename = settings.outputPrefix;
        filename += std::filesystem::path{"(" + tag + ").png"};
        film.save(filename);
    };

    const auto exportImage = [&](const std::string& tag, int shadeOptions) {
        film.postProcess(shadeOptions);
        saveCurrentImage(tag);
    };
    const auto exportFxaaPair = [&](const std::string& tag, int shadeOptions) {
        film.shade(shadeOptions);
        film.gammaCorrect();
        saveCurrentImage(tag);
        film.applyFxaa();
        saveCurrentImage(tag + "_FXAA");
    };
    const auto exportBeautyVariants =
        [&](const std::string& tag, int shadeOptions, bool depthOfField) {
        film.shade(shadeOptions);
        if (depthOfField) {
            film.depthOfFieldBlur();
        }
        const std::vector<vec3> linearPixels = film.pixels;

        film.gammaCorrect();
        saveCurrentImage(tag);
        film.applyFxaa();
        saveCurrentImage(tag + "_FXAA");

        film.pixels = linearPixels;
        film.bloom();
        film.gammaCorrect();
        saveCurrentImage(tag + "_Bloom");
        film.applyFxaa();
        saveCurrentImage(tag + "_Bloom_FXAA");
    };

    exportFxaaPair("DiffuseColor", Film::baseColor);
    exportImage("ShapeNormal", Film::shapeNormal);
    exportImage("SurfaceNormal", Film::surfaceNormal);

    film.spatialClamp();
    exportImage("Direct_Diffuse", Film::directDiffuse);
    exportImage("Direct_Specular", Film::directSpecular);
    exportImage("Indirect_Diffuse", Film::indirectDiffuse);
    exportImage("Indirect_Specular", Film::indirectSpecular);
    exportBeautyVariants("Raw", Film::full, false);

    film.filter();
    exportImage("Direct_Diffuse_Filter", Film::directDiffuse);
    exportImage("Direct_Specular_Filter", Film::directSpecular);
    exportImage("Indirect_Diffuse_Filter", Film::indirectDiffuse);
    exportImage("Indirect_Specular_Filter", Film::indirectSpecular);
    exportBeautyVariants("Filter", Film::full, false);

    if (settings.circleOfConfusion > kRayEpsilon) {
        film.focusDistance = settings.focusDistance;
        film.circleOfConfusion = settings.circleOfConfusion;
        film.cameraPosition = settings.position;
        exportImage("BaseColor_DepthOfField", Film::baseColor | Film::depthOfFieldEnabled);
        exportBeautyVariants("Filter_DepthOfField", Film::full, true);
    }
}

RenderStats renderToFiles(const Model& model, const RenderSettings& settings) {
    settings.validate();
    const auto totalStart = std::chrono::steady_clock::now();
    const int workerCount = settings.resolvedThreadCount();
    std::cout << "Rendering started with " << workerCount << " threads.\n";

    const std::filesystem::path parent = settings.outputPrefix.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    FilmRenderResult result = renderToFilmValidated(model, settings, workerCount);
    exportFilmToFiles(std::move(result.film), settings);

    const auto totalEnd = std::chrono::steady_clock::now();
    RenderStats stats = result.stats;
    stats.totalSeconds = std::chrono::duration<double>(totalEnd - totalStart).count();
    std::cout << "Rendering completed in " << stats.renderSeconds << " seconds.\n";
    std::cout << "Direct light samples: " << stats.directLightSamples << '\n';
    std::cout << "Post-processing finished. Total: " << stats.totalSeconds << " seconds.\n";
    return stats;
}

}  // namespace raym0nade
