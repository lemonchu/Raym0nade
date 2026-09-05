#include "raym0nade/sampling.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace raym0nade {

namespace sampling_detail {

std::size_t firstCumulativeWeightAbove(
    const std::vector<double>& cumulativeWeights, double target) noexcept {
    return static_cast<std::size_t>(
        std::upper_bound(
            cumulativeWeights.begin(), cumulativeWeights.end(), target) -
        cumulativeWeights.begin());
}

}  // namespace sampling_detail

namespace {

constexpr int kMaximumSamplingAttempts = 1;
constexpr float kThroughputLuminanceLimit = 64.0F;
constexpr float kMinimumPdf = 1.0e-8F;

float square(float value) noexcept {
    return value * value;
}

float fifthPower(float value) noexcept {
    const float squared = value * value;
    return squared * squared * value;
}

float schlickWeight(float cosine) noexcept {
    return fifthPower(std::clamp(1.0F - cosine, 0.0F, 1.0F));
}

float dielectricFresnel(float incidentCosine, float etaIncidentOverTransmitted) noexcept {
    if (!isFinite(incidentCosine) || !isFinite(etaIncidentOverTransmitted) ||
        etaIncidentOverTransmitted <= 0.0F) {
        return 1.0F;
    }

    incidentCosine = std::clamp(std::abs(incidentCosine), 0.0F, 1.0F);
    const float transmittedSineSquared = square(etaIncidentOverTransmitted) *
                                         std::max(0.0F, 1.0F - square(incidentCosine));
    if (transmittedSineSquared >= 1.0F) {
        return 1.0F;
    }

    const float transmittedCosine = safeSqrt(1.0F - transmittedSineSquared);
    const float perpendicularDenominator =
        etaIncidentOverTransmitted * incidentCosine + transmittedCosine;
    const float parallelDenominator =
        incidentCosine + etaIncidentOverTransmitted * transmittedCosine;
    if (perpendicularDenominator <= kMinimumPdf || parallelDenominator <= kMinimumPdf) {
        return 1.0F;
    }

    const float perpendicular =
        (etaIncidentOverTransmitted * incidentCosine - transmittedCosine) /
        perpendicularDenominator;
    const float parallel =
        (incidentCosine - etaIncidentOverTransmitted * transmittedCosine) /
        parallelDenominator;
    return std::clamp(0.5F * (square(perpendicular) + square(parallel)), 0.0F, 1.0F);
}

float gtr1(float normalDotHalf, float alpha) noexcept {
    alpha = std::clamp(alpha, 1.0e-3F, 0.9999F);
    const float alphaSquared = alpha * alpha;
    const float denominator = kPi * std::log(alphaSquared) *
                              (1.0F + (alphaSquared - 1.0F) * normalDotHalf * normalDotHalf);
    if (std::abs(denominator) <= kMinimumPdf) {
        return 0.0F;
    }
    return (alphaSquared - 1.0F) / denominator;
}

float gtr2(float normalDotHalf, float alpha) noexcept {
    alpha = std::max(alpha, 1.0e-3F);
    const float alphaSquared = alpha * alpha;
    const float term = 1.0F + (alphaSquared - 1.0F) * normalDotHalf * normalDotHalf;
    return alphaSquared / std::max(kPi * term * term, kMinimumPdf);
}

float smithGgx(float normalDotDirection, float alpha) noexcept {
    if (normalDotDirection <= 0.0F) {
        return 0.0F;
    }
    const float alphaSquared = alpha * alpha;
    const float cosineSquared = normalDotDirection * normalDotDirection;
    const float root = safeSqrt(alphaSquared + cosineSquared - alphaSquared * cosineSquared);
    return 1.0F / std::max(normalDotDirection + root, kMinimumPdf);
}

float cosinePdf(const vec3& normal, const vec3& direction) noexcept {
    return std::max(0.0F, glm::dot(normal, direction)) / kPi;
}

float microfacetReflectionPdf(
    const vec3& normal, const vec3& incoming, const vec3& outgoing, float roughness) noexcept {
    const vec3 halfVector = safeNormalize(incoming + outgoing);
    if (glm::dot(halfVector, halfVector) == 0.0F) {
        return 0.0F;
    }
    const float normalDotHalf = std::max(0.0F, glm::dot(normal, halfVector));
    const float outgoingDotHalf = std::abs(glm::dot(outgoing, halfVector));
    if (normalDotHalf <= 0.0F || outgoingDotHalf <= kMinimumPdf) {
        return 0.0F;
    }
    return gtr2(normalDotHalf, roughness) * normalDotHalf / (4.0F * outgoingDotHalf);
}

void limitThroughput(vec3& throughput) noexcept {
    const float luminance = glm::dot(throughput, kLuminanceWeights);
    if (isFinite(luminance) && luminance > kThroughputLuminanceLimit) {
        throughput *= kThroughputLuminanceLimit / luminance;
    }
}

vec3 refractDirection(const vec3& incoming, const vec3& normal, float eta) noexcept {
    const float normalDotIncoming = glm::dot(normal, incoming);
    const float discriminant = 1.0F - eta * eta * (1.0F - normalDotIncoming * normalDotIncoming);
    if (discriminant < 0.0F || !isFinite(discriminant)) {
        return vec3{std::numeric_limits<float>::quiet_NaN()};
    }
    return (eta * normalDotIncoming - safeSqrt(discriminant)) * normal - eta * incoming;
}

float triangleArea(const Face& face) noexcept {
    return 0.5F * glm::length(glm::cross(
                      face.vertices[1] - face.vertices[0], face.vertices[2] - face.vertices[0]));
}

vec3 sampleTriangle(
    const Face& face, Generator& generator, vec3& barycentricCoordinates) noexcept {
    const float root = safeSqrt(generator());
    const float second = generator();
    barycentricCoordinates = vec3{1.0F - root, root * (1.0F - second), root * second};
    return barycentricCoordinates[0] * face.vertices[0] +
           barycentricCoordinates[1] * face.vertices[1] +
           barycentricCoordinates[2] * face.vertices[2];
}

vec3 emissiveRadiance(const Face& face, const vec3& coordinates) {
    if (face.material == nullptr) {
        return vec3{0.0F};
    }
    const vec2 uv = coordinates[0] * face.vertexData[0]->uv +
                    coordinates[1] * face.vertexData[1]->uv +
                    coordinates[2] * face.vertexData[2]->uv;
    return face.material->emissiveColor(
        uv.x, uv.y, std::numeric_limits<float>::quiet_NaN());
}

double lightObjectWeight(const vec3& position, const LightObject& light) noexcept {
    const vec3 offset = light.center - position;
    const float distanceSquared = glm::dot(offset, offset);
    if (!isFinite(distanceSquared) || !isFinite(light.power) || light.power <= 0.0F) {
        return 0.0;
    }
    return static_cast<double>(light.power) /
           std::max(static_cast<double>(distanceSquared),
                    static_cast<double>(kMinimumLightDistanceSquared));
}

int sampleLightObject(
    const vec3& position,
    const Model& model,
    Generator& generator,
    float& objectProbability) noexcept {
    objectProbability = 0.0F;
    double totalWeight = 0.0;
    for (const LightObject& light : model.lights()) {
        totalWeight += lightObjectWeight(position, light);
    }
    if (!(totalWeight > 0.0) || !std::isfinite(totalWeight)) {
        return -1;
    }

    const double target = totalWeight * static_cast<double>(generator());
    double cumulative = 0.0;
    int lastValid = -1;
    double selectedWeight = 0.0;
    for (std::size_t index = 0; index < model.lights().size(); ++index) {
        const double weight = lightObjectWeight(position, model.lights()[index]);
        if (!(weight > 0.0)) {
            continue;
        }
        lastValid = static_cast<int>(index);
        cumulative += weight;
        if (target < cumulative) {
            selectedWeight = weight;
            objectProbability = static_cast<float>(selectedWeight / totalWeight);
            return static_cast<int>(index);
        }
    }

    if (lastValid >= 0) {
        selectedWeight = lightObjectWeight(
            position, model.lights()[static_cast<std::size_t>(lastValid)]);
        objectProbability = static_cast<float>(selectedWeight / totalWeight);
    }
    return lastValid;
}

struct AreaLightWeightCache {
    std::vector<double>& weights;
    std::vector<double>& cumulativeWeights;
    std::uint64_t& modelIdentity;
    vec3& position;
    std::size_t& lightCount;
    double& totalWeight;
    int& lastValidIndex;
    bool& initialized;
};

bool exactlyEqualPosition(const vec3& left, const vec3& right) noexcept {
    return left.x == right.x && left.y == right.y && left.z == right.z;
}

int sampleLightObjectCached(
    const vec3& position,
    const Model& model,
    Generator& generator,
    AreaLightWeightCache& cache,
    float& objectProbability) {
    const std::vector<LightObject>& lights = model.lights();
    if (!cache.initialized ||
        cache.modelIdentity != model.instanceIdentity() ||
        cache.lightCount != lights.size() ||
        !exactlyEqualPosition(cache.position, position)) {
        cache.weights.resize(lights.size());
        cache.cumulativeWeights.resize(lights.size());
        double totalWeight = 0.0;
        int lastValidIndex = -1;
        for (std::size_t index = 0; index < lights.size(); ++index) {
            const double weight = lightObjectWeight(position, lights[index]);
            cache.weights[index] = weight;
            totalWeight += weight;
            cache.cumulativeWeights[index] = totalWeight;
            if (weight > 0.0) {
                lastValidIndex = static_cast<int>(index);
            }
        }
        cache.modelIdentity = model.instanceIdentity();
        cache.position = position;
        cache.lightCount = lights.size();
        cache.totalWeight = totalWeight;
        cache.lastValidIndex = lastValidIndex;
        cache.initialized = true;
    }

    objectProbability = 0.0F;
    if (!(cache.totalWeight > 0.0) ||
        !std::isfinite(cache.totalWeight)) {
        return -1;
    }

    const double target =
        cache.totalWeight * static_cast<double>(generator());
    std::size_t selectedIndex = sampling_detail::firstCumulativeWeightAbove(
        cache.cumulativeWeights, target);
    if (selectedIndex >= cache.weights.size()) {
        if (cache.lastValidIndex < 0) {
            return -1;
        }
        selectedIndex =
            static_cast<std::size_t>(cache.lastValidIndex);
    }
    const double selectedWeight = cache.weights[selectedIndex];
    objectProbability =
        static_cast<float>(selectedWeight / cache.totalWeight);
    return static_cast<int>(selectedIndex);
}

bool sampleAreaLight(
    const Bsdf& bsdf,
    const Model& model,
    Generator& generator,
    float categoryProbability,
    AreaLightWeightCache* areaLightWeightCache,
    LightSample& result) {
    float objectProbability = 0.0F;
    const int objectIndex =
        areaLightWeightCache == nullptr
            ? sampleLightObject(
                  bsdf.surface.position, model, generator, objectProbability)
            : sampleLightObjectCached(
                  bsdf.surface.position,
                  model,
                  generator,
                  *areaLightWeightCache,
                  objectProbability);
    if (objectIndex < 0) {
        return false;
    }

    const LightObject& light = model.lights()[static_cast<std::size_t>(objectIndex)];
    const int faceIndex = light.faceDistribution.sample(generator);
    if (faceIndex < 0 || static_cast<std::size_t>(faceIndex) >= light.faces.size()) {
        return false;
    }
    const Face& face = light.faces[static_cast<std::size_t>(faceIndex)];
    const float area = triangleArea(face);
    const float faceProbability = light.faceDistribution.pdf(static_cast<std::size_t>(faceIndex));
    if (area <= kMinimumPdf || objectProbability <= 0.0F || faceProbability <= 0.0F) {
        return false;
    }

    vec3 coordinates;
    const vec3 lightPosition = sampleTriangle(face, generator, coordinates);
    const vec3 offset = lightPosition - bsdf.surface.position;
    const float distanceSquared = glm::dot(offset, offset);
    if (distanceSquared <= kMinimumLightDistanceSquared || !isFinite(distanceSquared)) {
        return false;
    }
    const float distance = safeSqrt(distanceSquared);
    const vec3 direction = offset / distance;
    const vec3 lightNormal = safeNormalize(glm::cross(
        face.vertices[1] - face.vertices[0], face.vertices[2] - face.vertices[0]));
    const float lightCosine = std::max(0.0F, glm::dot(lightNormal, -direction));
    if (lightCosine <= kMinimumPdf || model.occluded(
                                           Ray{bsdf.surface.position, direction},
                                           distance - kRayEpsilon)) {
        return false;
    }

    const float areaPdf = objectProbability * faceProbability / area;
    const float solidAnglePdf = areaPdf * distanceSquared / lightCosine;
    const float combinedPdf = categoryProbability * solidAnglePdf;
    if (combinedPdf <= kMinimumPdf || !isFinite(combinedPdf)) {
        return false;
    }

    result.throughput = bsdf.evaluate(direction);
    result.radiance = emissiveRadiance(face, coordinates) / combinedPdf;
    return isFinite(result.throughput) && isFinite(result.radiance);
}

bool sampleEnvironment(
    const Bsdf& bsdf,
    const Model& model,
    Generator& generator,
    float categoryProbability,
    LightSample& result) noexcept {
    int pixelIndex = -1;
    vec3 direction;
    vec3 weightedRadiance;
    if (!model.sky().sample(generator, pixelIndex, direction, weightedRadiance)) {
        return false;
    }
    if (model.occluded(
            Ray{bsdf.surface.position, direction}, std::numeric_limits<float>::infinity())) {
        return false;
    }
    result.throughput = bsdf.evaluate(direction);
    result.radiance = weightedRadiance / categoryProbability;
    return isFinite(result.throughput) && isFinite(result.radiance);
}

}  // namespace

Bsdf::Bsdf(const vec3& incoming, const vec3& hitPosition) : incomingDirection(incoming) {
    surface.position = hitPosition;
}

Bsdf::Bsdf(const vec3& incoming, const HitInfo& hitInfo)
    : incomingDirection(incoming), surface(hitInfo) {}

vec3 Bsdf::evaluateReflection(const vec3& outgoingDirection) const {
    const vec3& incoming = incomingDirection;
    const vec3& normal = surface.surfaceNormal;
    const float normalDotOutgoing = glm::dot(normal, outgoingDirection);
    const float normalDotIncoming = glm::dot(normal, incoming);
    if (normalDotOutgoing <= 0.0F || normalDotIncoming <= 0.0F) {
        return vec3{0.0F};
    }

    const vec3 halfVector = safeNormalize(outgoingDirection + incoming);
    if (glm::dot(halfVector, halfVector) == 0.0F) {
        return vec3{0.0F};
    }
    const float normalDotHalf = std::max(0.0F, glm::dot(normal, halfVector));
    const float outgoingDotHalf = std::clamp(glm::dot(outgoingDirection, halfVector), 0.0F, 1.0F);

    const vec3 baseColor = glm::max(surface.baseColor, vec3{0.0F});
    const vec3 specularColor = glm::mix(vec3{surface.specular}, baseColor, surface.metallic);

    const float outgoingFresnel = schlickWeight(normalDotOutgoing);
    const float incomingFresnel = schlickWeight(normalDotIncoming);
    const float diffuseAtGrazing = 0.5F + 2.0F * outgoingDotHalf * outgoingDotHalf * surface.roughness;
    const float diffuse = glm::mix(1.0F, diffuseAtGrazing, outgoingFresnel) *
                          glm::mix(1.0F, diffuseAtGrazing, incomingFresnel);

    const bool transmissive = surface.opacity < kRayEpsilon;
    const float fresnel = transmissive
                              ? dielectricFresnel(outgoingDotHalf, surface.eta)
                              : schlickWeight(outgoingDotHalf);
    const vec3 specularFresnel = transmissive
                                     ? vec3{fresnel}
                                     : glm::mix(specularColor, vec3{1.0F}, fresnel);
    const float distribution = gtr2(normalDotHalf, surface.roughness);
    const float geometry = smithGgx(normalDotOutgoing, surface.roughness) *
                           smithGgx(normalDotIncoming, surface.roughness);
    const float clearcoatDistribution = gtr1(normalDotHalf, glm::mix(0.1F, 0.001F, 0.2F));
    const float clearcoatGeometry = smithGgx(normalDotOutgoing, 0.25F) *
                                    smithGgx(normalDotIncoming, 0.25F);
    const float clearcoatFresnel = glm::mix(0.04F, 1.0F, fresnel);

    vec3 value = (1.0F / kPi) * diffuse * baseColor * (1.0F - surface.metallic) * surface.opacity;
    value += specularFresnel * distribution * geometry;
    value += surface.opacity * 0.375F * clearcoatGeometry * clearcoatFresnel *
             clearcoatDistribution * vec3{1.0F};
    value *= normalDotOutgoing;
    return isFinite(value) ? value : vec3{0.0F};
}

vec3 Bsdf::evaluateTransmission(const vec3& outgoingDirection) const {
    const vec3& normal = surface.surfaceNormal;
    const float normalDotOutgoing = glm::dot(normal, outgoingDirection);
    const float normalDotIncoming = glm::dot(normal, incomingDirection);
    if (normalDotOutgoing >= 0.0F || normalDotIncoming <= 0.0F) {
        return vec3{0.0F};
    }

    vec3 halfVector = safeNormalize(outgoingDirection + surface.eta * incomingDirection);
    if (glm::dot(halfVector, halfVector) == 0.0F) {
        return vec3{0.0F};
    }
    if (glm::dot(normal, halfVector) < 0.0F) {
        halfVector = -halfVector;
    }
    const float normalDotHalf = glm::dot(normal, halfVector);
    const float outgoingDotHalf = glm::dot(outgoingDirection, halfVector);
    const float incomingDotHalf = glm::dot(incomingDirection, halfVector);
    if (normalDotHalf <= 0.0F || outgoingDotHalf * normalDotOutgoing < 0.0F ||
        incomingDotHalf * normalDotIncoming < 0.0F) {
        return vec3{0.0F};
    }

    const float denominator = surface.eta * incomingDotHalf + outgoingDotHalf;
    if (std::abs(denominator) <= kMinimumPdf) {
        return vec3{0.0F};
    }

    const float distribution = gtr2(normalDotHalf, surface.roughness);
    const float smithProduct = smithGgx(std::abs(normalDotOutgoing), surface.roughness) *
                               smithGgx(std::abs(normalDotIncoming), surface.roughness);
    const float maskingShadowing = 4.0F * std::abs(normalDotOutgoing * normalDotIncoming) *
                                   smithProduct;
    const float fresnel = dielectricFresnel(incomingDotHalf, surface.eta);
    const float etaSquared = square(surface.eta);
    float value = distribution * (1.0F - fresnel) * maskingShadowing * etaSquared;
    value *= std::abs(outgoingDotHalf * incomingDotHalf) /
             std::max(std::abs(normalDotOutgoing * normalDotIncoming), kMinimumPdf);
    value *= -normalDotOutgoing / (denominator * denominator);
    return isFinite(value) ? vec3{value} : vec3{0.0F};
}

vec3 Bsdf::evaluate(const vec3& outgoingDirection) const {
    if (surface.opacity > 1.0F - kRayEpsilon || surface.entering) {
        return evaluateReflection(outgoingDirection);
    }
    return evaluateTransmission(outgoingDirection);
}

void Bsdf::sampleMicrofacet(
    Generator& generator, vec3& outgoingDirection, vec3& throughput, float& pdf, int& failures) const {
    vec3 tangent;
    vec3 bitangent;
    makeTangentSpace(surface.surfaceNormal, incomingDirection, tangent, bitangent);
    const float uniform = generator();
    const float azimuth = generator() * 2.0F * kPi;
    const float roughnessSquared = square(std::max(surface.roughness, 1.0e-3F));
    const float cosine = safeSqrt(
        (1.0F - uniform) / (1.0F + (roughnessSquared - 1.0F) * uniform));
    const float sine = safeSqrt(1.0F - cosine * cosine);
    const vec3 halfVector = sine * std::cos(azimuth) * tangent +
                            sine * std::sin(azimuth) * bitangent + cosine * surface.surfaceNormal;
    outgoingDirection = 2.0F * glm::dot(incomingDirection, halfVector) * halfVector - incomingDirection;
    pdf = microfacetReflectionPdf(
        surface.surfaceNormal, incomingDirection, outgoingDirection, surface.roughness);
    if (glm::dot(outgoingDirection, surface.shapeNormal) <= 0.0F || pdf <= 0.0F ||
        !isFinite(outgoingDirection)) {
        outgoingDirection = vec3{std::numeric_limits<float>::quiet_NaN()};
        throughput = vec3{0.0F};
        pdf = 0.0F;
        ++failures;
    }
}

void Bsdf::sampleCosine(
    Generator& generator, vec3& outgoingDirection, vec3& throughput, float& pdf, int& failures) const {
    vec3 tangent;
    vec3 bitangent;
    makeTangentSpace(surface.surfaceNormal, incomingDirection, tangent, bitangent);
    const float radius = safeSqrt(generator());
    const float azimuth = generator() * 2.0F * kPi;
    const float z = safeSqrt(1.0F - radius * radius);
    outgoingDirection = radius * std::cos(azimuth) * tangent + radius * std::sin(azimuth) * bitangent +
                        z * surface.surfaceNormal;
    pdf = cosinePdf(surface.surfaceNormal, outgoingDirection);
    if (glm::dot(outgoingDirection, surface.shapeNormal) <= 0.0F || pdf <= 0.0F ||
        !isFinite(outgoingDirection)) {
        outgoingDirection = vec3{std::numeric_limits<float>::quiet_NaN()};
        throughput = vec3{0.0F};
        pdf = 0.0F;
        ++failures;
    }
}

void Bsdf::sampleReflection(
    Generator& generator, vec3& direction, vec3& throughput, int& failures) const {
    constexpr float cosineTechniqueProbability = 0.5F;
    float selectedPdf = 0.0F;
    if (generator() < cosineTechniqueProbability) {
        sampleCosine(generator, direction, throughput, selectedPdf, failures);
    } else {
        sampleMicrofacet(generator, direction, throughput, selectedPdf, failures);
    }
    if (!isFinite(direction) || selectedPdf <= 0.0F) {
        throughput = vec3{0.0F};
        return;
    }

    const float diffusePdf = cosinePdf(surface.surfaceNormal, direction);
    const float specularPdf = microfacetReflectionPdf(
        surface.surfaceNormal, incomingDirection, direction, surface.roughness);
    const float mixturePdf = cosineTechniqueProbability * diffusePdf +
                             (1.0F - cosineTechniqueProbability) * specularPdf;
    if (mixturePdf <= kMinimumPdf || !isFinite(mixturePdf)) {
        direction = vec3{std::numeric_limits<float>::quiet_NaN()};
        throughput = vec3{0.0F};
        ++failures;
        return;
    }
    throughput = evaluateReflection(direction) / mixturePdf;
    limitThroughput(throughput);
}

void Bsdf::refractPrecisely(vec3& outgoingDirection, float& fresnel) const noexcept {
    const float normalDotIncoming = glm::dot(surface.surfaceNormal, incomingDirection);
    if (normalDotIncoming < 0.0F) {
        fresnel = 1.0F;
        outgoingDirection = vec3{std::numeric_limits<float>::quiet_NaN()};
        return;
    }
    outgoingDirection = refractDirection(incomingDirection, surface.surfaceNormal, surface.eta);
    if (!isFinite(outgoingDirection)) {
        fresnel = 1.0F;
        return;
    }
    fresnel = dielectricFresnel(normalDotIncoming, surface.eta);
}

void Bsdf::sampleTransmission(
    Generator& generator, vec3& direction, vec3& throughput, int& failures) const {
    if (std::abs(surface.eta - 1.0F) < kRayEpsilon) {
        direction = -incomingDirection;
        throughput = vec3{1.0F};
        return;
    }

    for (int attempt = 0; attempt < kMaximumSamplingAttempts; ++attempt) {
        vec3 tangent;
        vec3 bitangent;
        makeTangentSpace(surface.surfaceNormal, incomingDirection, tangent, bitangent);
        const float uniform = generator();
        const float azimuth = generator() * 2.0F * kPi;
        const float roughnessSquared = square(std::max(surface.roughness, 1.0e-3F));
        const float cosine = safeSqrt(
            (1.0F - uniform) / (1.0F + (roughnessSquared - 1.0F) * uniform));
        const float sine = safeSqrt(1.0F - cosine * cosine);
        const vec3 halfVector = sine * std::cos(azimuth) * tangent +
                                sine * std::sin(azimuth) * bitangent + cosine * surface.surfaceNormal;
        direction = refractDirection(incomingDirection, halfVector, surface.eta);
        if (!isFinite(direction) || glm::dot(direction, surface.shapeNormal) >= 0.0F) {
            ++failures;
            continue;
        }

        const float normalDotHalf = std::max(0.0F, glm::dot(surface.surfaceNormal, halfVector));
        const float incomingDotHalf = glm::dot(incomingDirection, halfVector);
        const float outgoingDotHalf = glm::dot(direction, halfVector);
        const float denominator = surface.eta * incomingDotHalf + outgoingDotHalf;
        if (std::abs(denominator) <= kMinimumPdf) {
            ++failures;
            continue;
        }
        const float halfPdf = gtr2(normalDotHalf, surface.roughness) * normalDotHalf;
        const float jacobian = std::abs(outgoingDotHalf) / (denominator * denominator);
        const float pdf = halfPdf * jacobian;
        if (pdf <= kMinimumPdf || !isFinite(pdf)) {
            ++failures;
            continue;
        }
        throughput = evaluateTransmission(direction) / pdf;
        limitThroughput(throughput);
        return;
    }

    direction = vec3{std::numeric_limits<float>::quiet_NaN()};
    throughput = vec3{0.0F};
}

namespace {

void sampleDirectLightImpl(
    const Bsdf& bsdf,
    const Model& model,
    Generator& generator,
    int sampleCount,
    std::vector<LightSample>& samples,
    AreaLightWeightCache* areaLightWeightCache) {
    samples.clear();
    if (sampleCount <= 0) {
        return;
    }
    samples.reserve(static_cast<std::size_t>(sampleCount));

    const bool hasEnvironment = !model.sky().empty();
    const bool hasAreaLights = !model.lights().empty();
    if (!hasEnvironment && !hasAreaLights) {
        return;
    }
    const float environmentProbability = hasEnvironment && hasAreaLights ? 0.5F : (hasEnvironment ? 1.0F : 0.0F);
    const float areaProbability = 1.0F - environmentProbability;

    for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex) {
        LightSample sample;
        sample.weight = 1.0F / static_cast<float>(sampleCount);
        const bool chooseEnvironment = hasEnvironment &&
                                       (!hasAreaLights || generator() < environmentProbability);
        const bool valid = chooseEnvironment
                               ? sampleEnvironment(bsdf, model, generator, environmentProbability, sample)
                               : sampleAreaLight(
                                     bsdf,
                                     model,
                                     generator,
                                     areaProbability,
                                     areaLightWeightCache,
                                     sample);
        if (valid) {
            limitThroughput(sample.throughput);
            samples.push_back(sample);
        }
    }
}

}  // namespace

void sampleDirectLight(
    const Bsdf& bsdf,
    const Model& model,
    Generator& generator,
    int sampleCount,
    std::vector<LightSample>& samples) {
    sampleDirectLightImpl(
        bsdf, model, generator, sampleCount, samples, nullptr);
}

void sampleDirectLight(
    const Bsdf& bsdf,
    const Model& model,
    Generator& generator,
    int sampleCount,
    std::vector<LightSample>& samples,
    DirectLightSamplingScratch& scratch) {
    AreaLightWeightCache cache{
        scratch.areaLightWeights_,
        scratch.areaLightCumulativeWeights_,
        scratch.areaLightModelIdentity_,
        scratch.areaLightPosition_,
        scratch.areaLightCount_,
        scratch.areaLightTotalWeight_,
        scratch.areaLightLastValidIndex_,
        scratch.areaLightWeightsInitialized_};
    sampleDirectLightImpl(
        bsdf,
        model,
        generator,
        sampleCount,
        samples,
        &cache);
}

}  // namespace raym0nade
