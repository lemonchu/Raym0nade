#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "raym0nade/gpu/vulkan_path_renderer.hpp"
#include "raym0nade/model.hpp"
#include "raym0nade/render.hpp"
#include "raym0nade/scene_data.hpp"

#ifndef RAYM0NADE_PATH_TRACE_SHADER
#error "RAYM0NADE_PATH_TRACE_SHADER must be defined by CMake."
#endif

#ifndef RAYM0NADE_TEST_SOURCE_DIR
#error "RAYM0NADE_TEST_SOURCE_DIR must be defined by CMake."
#endif

namespace {

using namespace raym0nade;
using namespace raym0nade::gpu;

constexpr std::string_view kNoCompatibleDevice =
    "No AMD Vulkan device satisfies the Ray Query backend requirements.";
constexpr float kComparisonAbsoluteTolerance = 2.0e-5F;
constexpr float kComparisonRelativeTolerance = 5.0e-5F;
constexpr vec3 kEnvironmentRadiance{0.25F, 0.5F, 1.0F};
constexpr vec3 kEmitterRadiance{8.0F, 6.0F, 4.0F};

void expect(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool exactlyEqual(const vec3& left, const vec3& right) noexcept {
    return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool near(float left, float right) noexcept {
    if (!std::isfinite(left) || !std::isfinite(right)) {
        return false;
    }
    const float scale = std::max(std::abs(left), std::abs(right));
    return std::abs(left - right) <=
           kComparisonAbsoluteTolerance + kComparisonRelativeTolerance * scale;
}

void expectNear(const vec3& actual, const vec3& expected, const std::string& message) {
    expect(
        isFinite(actual) && isFinite(expected) && near(actual.x, expected.x) &&
            near(actual.y, expected.y) && near(actual.z, expected.z),
        message);
}

float sumComponents(const vec3& value) noexcept {
    return value.x + value.y + value.z;
}

constexpr std::uint32_t packRgba8(
    std::uint32_t red,
    std::uint32_t green,
    std::uint32_t blue,
    std::uint32_t alpha) noexcept {
    return red | (green << 8U) | (blue << 16U) | (alpha << 24U);
}

PackedVertex makeVertex(
    const vec3& position, const vec3& normal, const vec2& uv) noexcept {
    return PackedVertex{
        {position.x, position.y, position.z, normal.x},
        {normal.y, normal.z, uv.x, uv.y},
    };
}

PackedMaterial makeMaterial(const vec3& diffuse, const vec3& emission) noexcept {
    PackedMaterial material;
    material.diffuseAndOpacity = {diffuse.x, diffuse.y, diffuse.z, 1.0F};
    material.emissionAndIor = {emission.x, emission.y, emission.z, 1.5F};
    material.transmissionAndRoughness = {1.0F, 1.0F, 1.0F, 0.5F};
    material.metallicSpecularAndReserved = {0.0F, 0.04F, 0.0F, 0.0F};
    return material;
}

std::array<std::uint32_t, 3> appendTriangle(
    PackedSceneData& scene,
    const std::array<vec3, 3>& positions,
    const vec3& normal,
    std::uint32_t materialId) {
    const std::uint32_t first = static_cast<std::uint32_t>(scene.vertices.size());
    constexpr std::array<vec2, 3> kUvs{
        vec2{0.0F, 0.0F},
        vec2{1.0F, 0.0F},
        vec2{0.5F, 1.0F},
    };
    for (std::size_t corner = 0U; corner < positions.size(); ++corner) {
        scene.vertices.push_back(makeVertex(positions[corner], normal, kUvs[corner]));
    }
    scene.triangleIndices.insert(
        scene.triangleIndices.end(), {first, first + 1U, first + 2U});
    scene.triangleMaterialIds.push_back(materialId);
    return {first, first + 1U, first + 2U};
}

std::uint32_t appendOnePixelTexture(
    PackedSceneData& scene, std::uint32_t texel) {
    const std::uint32_t textureId =
        static_cast<std::uint32_t>(scene.textures.size());
    const std::uint32_t mipId =
        static_cast<std::uint32_t>(scene.textureMipLevels.size());
    const std::uint32_t texelId =
        static_cast<std::uint32_t>(scene.textureTexelsRgba8.size());
    scene.textures.push_back(PackedTexture{mipId, 1U, 1U, 1U});
    scene.textureMipLevels.push_back(PackedTextureMip{texelId, 1U, 1U, 1U});
    scene.textureTexelsRgba8.push_back(texel);
    return textureId;
}

void addConstantEnvironment(PackedSceneData& scene) {
    scene.environment =
        PackedEnvironment{1U, 1U, kPackedEnvironmentHasImportance, 0U};
    scene.environmentRows = {
        PackedEnvironmentRow{1.0F, 1.0F, 4.0F * kPi, 0.0F},
    };
    scene.environmentTexels = {
        PackedEnvironmentTexel{{
            kEnvironmentRadiance.x,
            kEnvironmentRadiance.y,
            kEnvironmentRadiance.z,
            1.0F,
        }},
    };
}

PackedSceneData makeAreaLightScene(bool withEnvironment) {
    PackedSceneData scene;
    scene.materials.push_back(makeMaterial(vec3{0.7F, 0.4F, 0.2F}, vec3{0.0F}));
    scene.materials.push_back(makeMaterial(vec3{1.0F}, kEmitterRadiance));

    appendTriangle(
        scene,
        {
            vec3{-10.0F, -10.0F, 3.0F},
            vec3{10.0F, -10.0F, 3.0F},
            vec3{0.0F, 10.0F, 3.0F},
        },
        vec3{0.0F, 0.0F, -1.0F},
        0U);
    const std::array<std::uint32_t, 3> emitterVertices = appendTriangle(
        scene,
        {
            vec3{1.0F, -0.5F, 1.0F},
            vec3{2.0F, -0.5F, 1.0F},
            vec3{1.5F, 0.5F, 1.0F},
        },
        vec3{0.0F, 0.0F, 1.0F},
        1U);

    scene.areaLights.push_back(PackedAreaLight{
        {1.5F, -1.0F / 6.0F, 1.0F, 3.2F},
        {0U, 1U, 0U, 0U},
    });
    scene.areaLightTriangles.push_back(PackedAreaLightTriangle{
        {emitterVertices[0], emitterVertices[1], emitterVertices[2], 1U},
        {0.5F, 1.0F, 1.0F, 0.0F},
    });
    if (withEnvironment) {
        addConstantEnvironment(scene);
    }
    scene.validate();
    return scene;
}

PackedSceneData makeMaterialAndTransmissionScene() {
    PackedSceneData scene;

    PackedMaterial transparent = makeMaterial(vec3{0.9F}, vec3{0.0F});
    transparent.diffuseAndOpacity[3] = 0.0F;
    transparent.emissionAndIor[3] = 1.5F;
    transparent.transmissionAndRoughness = {0.9F, 0.9F, 0.9F, 0.001F};
    scene.materials.push_back(transparent);
    scene.materials.push_back(makeMaterial(vec3{0.65F, 0.7F, 0.75F}, vec3{0.0F}));
    scene.materials.back().transmissionAndRoughness[3] = 1.0F;
    scene.materials.push_back(makeMaterial(vec3{1.0F}, vec3{20.0F, 16.0F, 12.0F}));
    scene.materials.push_back(makeMaterial(vec3{1.0F}, vec3{0.0F}));
    scene.materials.push_back(
        makeMaterial(vec3{0.5F, 0.75F, 1.0F}, vec3{2.0F, 3.0F, 4.0F}));

    const std::uint32_t cutoutTexture = appendOnePixelTexture(
        scene, packRgba8(255U, 0U, 0U, 0U));
    const std::uint32_t diffuseTexture = appendOnePixelTexture(
        scene, packRgba8(128U, 64U, 32U, 255U));
    const std::uint32_t specularTexture = appendOnePixelTexture(
        scene, packRgba8(0U, 64U, 192U, 255U));
    const std::uint32_t emissiveTexture = appendOnePixelTexture(
        scene, packRgba8(128U, 128U, 128U, 255U));
    const std::uint32_t normalTexture = appendOnePixelTexture(
        scene, packRgba8(191U, 128U, 238U, 255U));

    PackedMaterial& cutout = scene.materials[3];
    cutout.flagsAndReserved[0] =
        kPackedMaterialHasDiffuseTexture | kPackedMaterialCutout;
    cutout.textureIds[0] = cutoutTexture;

    PackedMaterial& textured = scene.materials[4];
    textured.flagsAndReserved[0] =
        kPackedMaterialHasDiffuseTexture | kPackedMaterialHasSpecularTexture |
        kPackedMaterialHasEmissiveTexture | kPackedMaterialHasNormalTexture;
    textured.textureIds = {
        diffuseTexture,
        specularTexture,
        emissiveTexture,
        normalTexture,
    };

    appendTriangle(
        scene,
        {
            vec3{-0.5F, -0.5F, 1.0F},
            vec3{0.0F, 0.5F, 1.0F},
            vec3{0.5F, -0.5F, 1.0F},
        },
        vec3{0.0F, 0.0F, -1.0F},
        0U);
    appendTriangle(
        scene,
        {
            vec3{-0.5F, -0.5F, 1.2F},
            vec3{0.5F, -0.5F, 1.2F},
            vec3{0.0F, 0.5F, 1.2F},
        },
        vec3{0.0F, 0.0F, 1.0F},
        0U);
    appendTriangle(
        scene,
        {
            vec3{-10.0F, -10.0F, 2.5F},
            vec3{0.0F, 10.0F, 2.5F},
            vec3{10.0F, -10.0F, 2.5F},
        },
        vec3{0.0F, 0.0F, -1.0F},
        1U);
    const std::array<std::uint32_t, 3> emitterVertices = appendTriangle(
        scene,
        {
            vec3{1.5F, -0.5F, 0.5F},
            vec3{2.5F, -0.5F, 0.5F},
            vec3{2.0F, 0.5F, 0.5F},
        },
        vec3{0.0F, 0.0F, 1.0F},
        2U);
    appendTriangle(
        scene,
        {
            vec3{-3.5F, -0.5F, 1.0F},
            vec3{-3.0F, 0.5F, 1.0F},
            vec3{-2.5F, -0.5F, 1.0F},
        },
        vec3{0.0F, 0.0F, -1.0F},
        3U);
    appendTriangle(
        scene,
        {
            vec3{-3.75F, -0.75F, 1.5F},
            vec3{-3.0F, 0.75F, 1.5F},
            vec3{-2.25F, -0.75F, 1.5F},
        },
        vec3{0.0F, 0.0F, -1.0F},
        4U);

    scene.areaLights.push_back(PackedAreaLight{
        {2.0F, -1.0F / 6.0F, 0.5F, 8.0F},
        {0U, 1U, 0U, 0U},
    });
    scene.areaLightTriangles.push_back(PackedAreaLightTriangle{
        {emitterVertices[0], emitterVertices[1], emitterVertices[2], 2U},
        {0.5F, 1.0F, 1.0F, 0.0F},
    });
    scene.validate();
    return scene;
}

RenderSettings makeSettings(int width, int height, int samplesPerPixel) {
    RenderSettings settings;
    settings.width = width;
    settings.height = height;
    settings.samplesPerPixel = samplesPerPixel;
    settings.position = vec3{0.0F};
    settings.direction = vec3{0.0F, 0.0F, 1.0F};
    settings.up = vec3{0.0F, 1.0F, 0.0F};
    settings.right = vec3{1.0F, 0.0F, 0.0F};
    settings.pixelScale = 0.1F;
    settings.exposure = 1.0F;
    settings.directLightProbability = 1.0F;
    settings.seed = 0x5a17c3e9U;
    return settings;
}

void validateRadiance(const RadianceData& value, const std::string& description) {
    expect(isFinite(value.radiance), description + " radiance is not finite.");
    expect(
        std::isfinite(value.varianceAccumulator) &&
            value.varianceAccumulator >= 0.0F,
        description + " variance is invalid.");
}

void validateFilm(
    const Film& film, bool expectEveryPixelToHit, const std::string& description) {
    const std::size_t pixelCount =
        static_cast<std::size_t>(film.width()) * static_cast<std::size_t>(film.height());
    expect(film.width() > 0 && film.height() > 0, description + " extent is invalid.");
    expect(film.gBuffer.size() == pixelCount, description + " GBuffer size is invalid.");
    expect(
        film.directDiffuseRadiance.size() == pixelCount &&
            film.directSpecularRadiance.size() == pixelCount &&
            film.indirectDiffuseRadiance.size() == pixelCount &&
            film.indirectSpecularRadiance.size() == pixelCount &&
            film.pixels.size() == pixelCount,
        description + " Film plane sizes are invalid.");

    for (std::size_t index = 0U; index < pixelCount; ++index) {
        const HitInfo& hit = film.gBuffer[index];
        if (expectEveryPixelToHit) {
            expect(
                isFinite(hit.shapeNormal) && isFinite(hit.surfaceNormal) &&
                    isFinite(hit.position),
                description + " contains a non-finite hit geometry value.");
        } else if (!isFinite(hit.position)) {
            expect(
                !isFinite(hit.shapeNormal) && !isFinite(hit.surfaceNormal),
                description + " miss sentinels are inconsistent.");
        }
        expect(
            isFinite(hit.emission) && isFinite(hit.baseColor) &&
                std::isfinite(hit.specular) && std::isfinite(hit.roughness) &&
                std::isfinite(hit.metallic) && std::isfinite(hit.opacity) &&
                std::isfinite(hit.eta),
            description + " contains a non-finite material value.");
        validateRadiance(
            film.directDiffuseRadiance[index], description + " direct diffuse");
        validateRadiance(
            film.directSpecularRadiance[index], description + " direct specular");
        validateRadiance(
            film.indirectDiffuseRadiance[index], description + " indirect diffuse");
        validateRadiance(
            film.indirectSpecularRadiance[index], description + " indirect specular");
        expect(isFinite(film.pixels[index]), description + " display pixel is not finite.");
    }
}

void expectSameHit(const HitInfo& left, const HitInfo& right, const std::string& message) {
    expect(
        exactlyEqual(left.shapeNormal, right.shapeNormal) &&
            exactlyEqual(left.surfaceNormal, right.surfaceNormal) &&
            exactlyEqual(left.emission, right.emission) &&
            exactlyEqual(left.baseColor, right.baseColor) &&
            exactlyEqual(left.position, right.position) &&
            left.specular == right.specular && left.roughness == right.roughness &&
            left.metallic == right.metallic && left.opacity == right.opacity &&
            left.eta == right.eta && left.materialId == right.materialId &&
            left.entering == right.entering,
        message);
}

void expectSameRadiance(
    const RadianceData& left,
    const RadianceData& right,
    bool exact,
    const std::string& message) {
    if (exact) {
        expect(
            exactlyEqual(left.radiance, right.radiance) &&
                left.varianceAccumulator == right.varianceAccumulator,
            message);
        return;
    }
    expectNear(left.radiance, right.radiance, message + " radiance differs.");
    expect(
        near(left.varianceAccumulator, right.varianceAccumulator),
        message + " variance differs.");
}

void compareFilm(
    const Film& left,
    const Film& right,
    bool exactRadiance,
    const std::string& description) {
    expect(
        left.width() == right.width() && left.height() == right.height(),
        description + " extent differs.");
    expect(left.gBuffer.size() == right.gBuffer.size(), description + " size differs.");
    for (std::size_t index = 0U; index < left.gBuffer.size(); ++index) {
        expectSameHit(
            left.gBuffer[index],
            right.gBuffer[index],
            description + " GBuffer differs at pixel " + std::to_string(index) + '.');
        expectSameRadiance(
            left.directDiffuseRadiance[index],
            right.directDiffuseRadiance[index],
            exactRadiance,
            description + " direct diffuse differs at pixel " +
                std::to_string(index) + '.');
        expectSameRadiance(
            left.directSpecularRadiance[index],
            right.directSpecularRadiance[index],
            exactRadiance,
            description + " direct specular differs at pixel " +
                std::to_string(index) + '.');
        expectSameRadiance(
            left.indirectDiffuseRadiance[index],
            right.indirectDiffuseRadiance[index],
            exactRadiance,
            description + " indirect diffuse differs at pixel " +
                std::to_string(index) + '.');
        expectSameRadiance(
            left.indirectSpecularRadiance[index],
            right.indirectSpecularRadiance[index],
            exactRadiance,
            description + " indirect specular differs at pixel " +
                std::to_string(index) + '.');
        expect(
            exactlyEqual(left.pixels[index], right.pixels[index]),
            description + " display pixel differs at pixel " +
                std::to_string(index) + '.');
    }
}

void validateTimings(
    const VulkanPathRenderTimings& timings,
    std::uint64_t expectedDispatches,
    const std::string& description) {
    expect(
        std::isfinite(timings.hostRenderMilliseconds) &&
            timings.hostRenderMilliseconds >= 0.0,
        description + " host timing is invalid.");
    expect(
        timings.dispatchCount == expectedDispatches,
        description + " dispatch count is incorrect.");
    if (timings.gpuTimestampAvailable) {
        expect(
            std::isfinite(timings.gpuDispatchMilliseconds) &&
                timings.gpuDispatchMilliseconds >= 0.0,
            description + " GPU timing is invalid.");
    }
}

void validateStats(const RenderStats& stats, const std::string& description) {
    expect(
        std::isfinite(stats.renderSeconds) && stats.renderSeconds >= 0.0 &&
            std::isfinite(stats.totalSeconds) && stats.totalSeconds >= 0.0,
        description + " render statistics are invalid.");
}

void expectCleanValidation(
    const VulkanValidationReport& report, const std::string& description) {
    expect(
        report.requested && report.enabled &&
            report.synchronizationValidationEnabled && report.errorCount == 0U &&
            report.warningCount == 0U,
        description + " did not run clean validation and synchronization validation; it reported " +
            std::to_string(report.errorCount) +
            " validation error(s) and " + std::to_string(report.warningCount) +
            " warning(s).");
}

void expectStatisticallyNear(
    const vec3& actual,
    const vec3& expected,
    float relativeTolerance,
    const std::string& description) {
    expect(isFinite(actual) && isFinite(expected), description + " is not finite.");
    for (int component = 0; component < 3; ++component) {
        const float scale = std::max(std::abs(expected[component]), 1.0e-3F);
        expect(
            std::abs(actual[component] - expected[component]) <=
                relativeTolerance * scale + 2.0e-4F,
            description + " differs beyond its statistical tolerance in component " +
                std::to_string(component) + '.');
    }
}

VulkanPathRenderResult renderOnce(
    const PackedSceneData& scene,
    const std::filesystem::path& shaderPath,
    const RenderSettings& settings,
    VulkanPathRenderOptions options,
    std::uint64_t expectedDispatches,
    const std::string& description) {
    VulkanPathRenderer renderer{scene, shaderPath, options};
    VulkanPathRenderResult result = renderer.render(settings);
    validateTimings(result.timings, expectedDispatches, description);
    validateStats(result.stats, description);
    expectCleanValidation(renderer.validationReport(), description);
    return result;
}

void testDeterminismAndPartitionInvariance(
    const std::filesystem::path& shaderPath,
    const VulkanRayQueryOptions& vulkanOptions) {
    const PackedSceneData scene = makeAreaLightScene(false);
    const RenderSettings settings = makeSettings(5, 4, 7);
    constexpr std::uint64_t kExpectedDirectSamples = 5U * 4U * 7U;

    VulkanPathRenderOptions referenceOptions;
    referenceOptions.vulkan = vulkanOptions;
    referenceOptions.tileWidth = 16U;
    referenceOptions.tileHeight = 16U;
    referenceOptions.samplesPerBatch = 7U;
    VulkanPathRenderResult reference = [&] {
        VulkanPathRenderer renderer{scene, shaderPath, referenceOptions};
        VulkanPathRenderResult first = renderer.render(settings);
        VulkanPathRenderResult repeated = renderer.render(settings);
        validateFilm(first.film, true, "Reference path render");
        validateFilm(repeated.film, true, "Repeated path render");
        compareFilm(first.film, repeated.film, true, "Repeated path render");
        expect(
            first.stats.directLightSamples == kExpectedDirectSamples &&
                repeated.stats.directLightSamples == kExpectedDirectSamples,
            "A direct-only render must count one direct sample per pixel and SPP.");
        validateTimings(first.timings, 1U, "Reference path render");
        validateTimings(repeated.timings, 1U, "Repeated path render");
        validateStats(first.stats, "Reference path render");
        validateStats(repeated.stats, "Repeated path render");
        expectCleanValidation(renderer.validationReport(), "Reference path renderer");
        std::cout << "GPU path test device: " << renderer.deviceName() << '\n';
        return first;
    }();

    const std::size_t centerIndex = 2U * 5U + 2U;
    const vec3 centerDirect =
        reference.film.directDiffuseRadiance[centerIndex].radiance +
        reference.film.directSpecularRadiance[centerIndex].radiance;
    expect(
        isFinite(centerDirect) && sumComponents(centerDirect) > 1.0e-4F,
        "The unobstructed receiver must obtain positive direct area-light radiance.");

    VulkanPathRenderOptions tiledOptions = referenceOptions;
    tiledOptions.tileWidth = 2U;
    tiledOptions.tileHeight = 3U;
    VulkanPathRenderResult tiled = renderOnce(
        scene,
        shaderPath,
        settings,
        tiledOptions,
        6U,
        "Tiled path render");
    validateFilm(tiled.film, true, "Tiled path render");
    compareFilm(reference.film, tiled.film, true, "Tile partition invariance");
    expect(
        tiled.stats.directLightSamples == kExpectedDirectSamples,
        "Tile partitioning changed the direct-light sample count.");

    VulkanPathRenderOptions batchedOptions = referenceOptions;
    batchedOptions.samplesPerBatch = 2U;
    VulkanPathRenderResult batched = renderOnce(
        scene,
        shaderPath,
        settings,
        batchedOptions,
        4U,
        "Batched path render");
    validateFilm(batched.film, true, "Batched path render");
    compareFilm(reference.film, batched.film, true, "Sample-batch invariance");
    expect(
        batched.stats.directLightSamples == kExpectedDirectSamples,
        "Sample batching changed the direct-light sample count.");
    batched.film.shade(Film::full);
    expect(
        isFinite(batched.film.pixels[centerIndex]) &&
            sumComponents(batched.film.pixels[centerIndex]) > 1.0e-4F,
        "Film shading did not expose the positive area-light contribution.");
}

void testEnvironmentMissAndEmissiveHit(
    const std::filesystem::path& shaderPath,
    const VulkanRayQueryOptions& vulkanOptions) {
    const PackedSceneData scene = makeAreaLightScene(true);
    VulkanPathRenderOptions options;
    options.vulkan = vulkanOptions;
    options.tileWidth = 2U;
    options.tileHeight = 2U;
    options.samplesPerBatch = 2U;
    VulkanPathRenderer renderer{scene, shaderPath, options};

    RenderSettings indirectSettings = makeSettings(3, 2, 5);
    indirectSettings.directLightProbability = 0.0F;
    const VulkanPathRenderResult indirect = renderer.render(indirectSettings);
    validateFilm(indirect.film, true, "Indirect-path smoke render");
    validateTimings(indirect.timings, 6U, "Indirect-path smoke render");
    validateStats(indirect.stats, "Indirect-path smoke render");
    for (std::size_t index = 0U; index < indirect.film.gBuffer.size(); ++index) {
        expect(
            exactlyEqual(
                indirect.film.directDiffuseRadiance[index].radiance, vec3{0.0F}) &&
                exactlyEqual(
                    indirect.film.directSpecularRadiance[index].radiance, vec3{0.0F}),
            "An indirect-only render wrote first-hit direct radiance.");
    }

    RenderSettings missSettings = makeSettings(2, 2, 3);
    missSettings.position = vec3{20.0F, 0.0F, 0.0F};
    missSettings.directLightProbability = 0.5F;
    VulkanPathRenderResult miss = renderer.render(missSettings);
    validateFilm(miss.film, false, "Environment-miss path render");
    validateTimings(miss.timings, 2U, "Environment-miss path render");
    validateStats(miss.stats, "Environment-miss path render");
    expect(
        miss.stats.directLightSamples == 0U,
        "Primary misses must not count direct-light samples.");
    for (std::size_t index = 0U; index < miss.film.gBuffer.size(); ++index) {
        const HitInfo& hit = miss.film.gBuffer[index];
        expect(!isFinite(hit.position), "The environment fixture unexpectedly hit geometry.");
        expectNear(
            hit.emission,
            kEnvironmentRadiance,
            "A primary miss did not preserve constant environment radiance.");
    }
    miss.film.shade(Film::full);
    for (const vec3& pixel : miss.film.pixels) {
        expectNear(
            pixel,
            kEnvironmentRadiance,
            "Film emission shading did not expose environment radiance.");
    }

    RenderSettings emissionSettings = makeSettings(2, 2, 3);
    emissionSettings.position = vec3{1.5F, -1.0F / 6.0F, 0.0F};
    const VulkanPathRenderResult emission = renderer.render(emissionSettings);
    validateFilm(emission.film, true, "Emissive-hit path render");
    validateTimings(emission.timings, 2U, "Emissive-hit path render");
    validateStats(emission.stats, "Emissive-hit path render");
    constexpr std::size_t kCenterPixel = 3U;
    const HitInfo& center = emission.film.gBuffer[kCenterPixel];
    expect(center.materialId == 1, "The emissive center ray hit the wrong material.");
    expectNear(
        center.emission,
        kEmitterRadiance,
        "A primary emissive hit returned the wrong emission.");
    expect(
        exactlyEqual(emission.film.directDiffuseRadiance[kCenterPixel].radiance, vec3{0.0F}) &&
            exactlyEqual(emission.film.directSpecularRadiance[kCenterPixel].radiance, vec3{0.0F}) &&
            exactlyEqual(emission.film.indirectDiffuseRadiance[kCenterPixel].radiance, vec3{0.0F}) &&
            exactlyEqual(emission.film.indirectSpecularRadiance[kCenterPixel].radiance, vec3{0.0F}),
        "A directly visible emitter must be represented by GBuffer emission, not path radiance.");
    expectCleanValidation(renderer.validationReport(), "Environment/emission path renderer");
}

float decodeColor(std::uint32_t encoded) {
    return std::pow(static_cast<float>(encoded) / 255.0F, 2.2F);
}

void testMaterialTexturesCutoutAndTransmission(
    const std::filesystem::path& shaderPath,
    const VulkanRayQueryOptions& vulkanOptions) {
    const PackedSceneData scene = makeMaterialAndTransmissionScene();
    VulkanPathRenderOptions options;
    options.vulkan = vulkanOptions;
    options.vulkan.textureTexelPageBytes = 16U;
    options.tileWidth = 2U;
    options.tileHeight = 2U;
    options.samplesPerBatch = 8U;
    VulkanPathRenderer renderer{scene, shaderPath, options};

    RenderSettings textureSettings = makeSettings(2, 2, 1);
    textureSettings.position = vec3{-3.0F, 0.0F, 0.0F};
    const VulkanPathRenderResult textureResult = renderer.render(textureSettings);
    validateFilm(textureResult.film, true, "Textured-cutout path render");
    validateTimings(textureResult.timings, 1U, "Textured-cutout path render");
    validateStats(textureResult.stats, "Textured-cutout path render");

    const vec3 expectedBaseColor{
        0.5F * decodeColor(128U),
        0.75F * decodeColor(64U),
        decodeColor(32U),
    };
    const vec3 expectedEmission{
        2.0F * decodeColor(128U),
        3.0F * decodeColor(128U),
        4.0F * decodeColor(128U),
    };
    for (const HitInfo& hit : textureResult.film.gBuffer) {
        expect(
            hit.materialId == 4,
            "The zero-alpha front candidate did not reveal the textured rear material.");
        expectNear(hit.baseColor, expectedBaseColor, "Diffuse texture decoding differs.");
        expectNear(hit.emission, expectedEmission, "Emissive texture decoding differs.");
        expect(
            near(hit.roughness, 64.0F / 255.0F) &&
                near(hit.metallic, 192.0F / 255.0F),
            "The specular texture did not provide roughness and metallic values.");
        const vec3 normalDifference = hit.surfaceNormal - hit.shapeNormal;
        expect(
            isFinite(hit.surfaceNormal) &&
                glm::dot(normalDifference, normalDifference) > 1.0e-2F &&
                glm::dot(hit.surfaceNormal, hit.shapeNormal) > 0.0F,
            "The normal texture did not produce a finite, front-facing perturbation.");
    }

    RenderSettings transmissionSettings = makeSettings(2, 2, 64);
    transmissionSettings.pixelScale = 0.05F;
    transmissionSettings.directLightProbability = 0.0F;
    const VulkanPathRenderResult transmission = renderer.render(transmissionSettings);
    const VulkanPathRenderResult transmissionRepeated = renderer.render(transmissionSettings);
    validateFilm(transmission.film, true, "Transmission path render");
    validateFilm(transmissionRepeated.film, true, "Repeated transmission path render");
    validateTimings(transmission.timings, 8U, "Transmission path render");
    validateTimings(
        transmissionRepeated.timings, 8U, "Repeated transmission path render");
    validateStats(transmission.stats, "Transmission path render");
    validateStats(transmissionRepeated.stats, "Repeated transmission path render");
    compareFilm(
        transmission.film,
        transmissionRepeated.film,
        true,
        "Repeated transmission path render");

    float totalIndirectSpecular = 0.0F;
    for (std::size_t index = 0U; index < transmission.film.gBuffer.size(); ++index) {
        const HitInfo& hit = transmission.film.gBuffer[index];
        expect(
            hit.materialId == 0 && hit.entering && hit.opacity == 0.0F &&
                hit.eta == 1.5F,
            "The primary transparent boundary GBuffer contract changed.");
        expect(
            exactlyEqual(
                transmission.film.directDiffuseRadiance[index].radiance, vec3{0.0F}) &&
                exactlyEqual(
                    transmission.film.directSpecularRadiance[index].radiance, vec3{0.0F}),
            "An indirect-only transmission render wrote first-hit direct radiance.");
        totalIndirectSpecular += sumComponents(
            transmission.film.indirectSpecularRadiance[index].radiance);
    }
    expect(
        transmission.stats.directLightSamples > 0U,
        "The two-boundary indirect paths never reached terminal direct-light sampling.");
    expect(
        std::isfinite(totalIndirectSpecular) && totalIndirectSpecular > 1.0e-4F,
        "Transparent enter/exit paths did not carry finite area-light radiance.");

    PackedSceneData mixedScene = scene;
    mixedScene.materials[0].transmissionAndRoughness[3] = 0.5F;
    VulkanPathRenderer mixedRenderer{mixedScene, shaderPath, options};
    RenderSettings mixedSettings = transmissionSettings;
    mixedSettings.samplesPerPixel = 128;
    mixedSettings.directLightProbability = 0.35F;
    const VulkanPathRenderResult mixed = mixedRenderer.render(mixedSettings);
    validateFilm(mixed.film, true, "Mixed-technique transmission path render");
    validateTimings(mixed.timings, 16U, "Mixed-technique transmission path render");
    validateStats(mixed.stats, "Mixed-technique transmission path render");
    expect(
        mixed.stats.directLightSamples > 0U,
        "The mixed direct/indirect render did not record any direct-light decisions.");
    float mixedDirect = 0.0F;
    float mixedIndirect = 0.0F;
    for (std::size_t index = 0U; index < mixed.film.gBuffer.size(); ++index) {
        mixedDirect += sumComponents(
            mixed.film.directDiffuseRadiance[index].radiance +
            mixed.film.directSpecularRadiance[index].radiance);
        mixedIndirect += sumComponents(
            mixed.film.indirectDiffuseRadiance[index].radiance +
            mixed.film.indirectSpecularRadiance[index].radiance);
    }
    expect(
        std::isfinite(mixedDirect) && mixedDirect > 1.0e-4F,
        "The mixed-technique render did not preserve any direct-light radiance.");
    expect(
        std::isfinite(mixedIndirect) && mixedIndirect > 1.0e-4F,
        "The mixed-technique render did not preserve any indirect-light radiance.");
    expectCleanValidation(mixedRenderer.validationReport(), "Mixed transmission path renderer");
    expectCleanValidation(renderer.validationReport(), "Material/transmission path renderer");
}

void testImportedCpuReference(
    const std::filesystem::path& shaderPath,
    const VulkanRayQueryOptions& vulkanOptions) {
    const std::filesystem::path sourceDirectory{RAYM0NADE_TEST_SOURCE_DIR};
    const Model model{
        sourceDirectory / "tests" / "data", "triangle.obj", "null"};
    expect(
        model.faceCount() == 2U && model.lights().size() == 1U,
        "The imported CPU/GPU path fixture is incomplete.");

    RenderSettings settings = makeSettings(1, 1, 32768);
    settings.threadCount = 1;
    settings.directLightProbability = 1.0F;
    const FilmRenderResult cpu = renderToFilm(model, settings);
    validateFilm(cpu.film, true, "Imported CPU reference");
    expect(
        cpu.stats.directLightSamples ==
            static_cast<std::uint64_t>(settings.samplesPerPixel),
        "The imported CPU reference did not execute one direct sample per SPP.");

    VulkanPathRenderOptions options;
    options.vulkan = vulkanOptions;
    options.tileWidth = 1U;
    options.tileHeight = 1U;
    options.samplesPerBatch = 64U;
    VulkanPathRenderer renderer{model.packScene(), shaderPath, options};
    const VulkanPathRenderResult gpu = renderer.render(settings);
    validateFilm(gpu.film, true, "Imported GPU comparison");
    validateTimings(gpu.timings, 512U, "Imported GPU comparison");
    validateStats(gpu.stats, "Imported GPU comparison");
    expect(
        gpu.stats.directLightSamples == cpu.stats.directLightSamples,
        "The imported CPU/GPU comparison executed different direct sample counts.");

    const HitInfo& cpuHit = cpu.film.gBuffer.front();
    const HitInfo& gpuHit = gpu.film.gBuffer.front();
    expect(
        cpuHit.materialId == gpuHit.materialId && cpuHit.entering == gpuHit.entering,
        "The imported CPU/GPU comparison classified the primary hit differently.");
    expectNear(gpuHit.position, cpuHit.position, "Imported CPU/GPU hit position differs.");
    expectNear(gpuHit.shapeNormal, cpuHit.shapeNormal, "Imported CPU/GPU shape normal differs.");
    expectNear(
        gpuHit.surfaceNormal,
        cpuHit.surfaceNormal,
        "Imported CPU/GPU surface normal differs.");
    expectNear(gpuHit.baseColor, cpuHit.baseColor, "Imported CPU/GPU base color differs.");
    expectNear(gpuHit.emission, cpuHit.emission, "Imported CPU/GPU emission differs.");

    const vec3 cpuDirectDiffuse =
        cpu.film.directDiffuseRadiance.front().radiance;
    const vec3 cpuDirectSpecular =
        cpu.film.directSpecularRadiance.front().radiance;
    const vec3 gpuDirectDiffuse =
        gpu.film.directDiffuseRadiance.front().radiance;
    const vec3 gpuDirectSpecular =
        gpu.film.directSpecularRadiance.front().radiance;
    expectStatisticallyNear(
        gpuDirectDiffuse,
        cpuDirectDiffuse,
        0.06F,
        "Imported CPU/GPU direct-diffuse estimator");
    expectStatisticallyNear(
        gpuDirectSpecular,
        cpuDirectSpecular,
        0.12F,
        "Imported CPU/GPU direct-specular estimator");
    expectStatisticallyNear(
        gpuDirectDiffuse + gpuDirectSpecular,
        cpuDirectDiffuse + cpuDirectSpecular,
        0.05F,
        "Imported CPU/GPU total direct estimator");
    expectNear(
        gpu.film.indirectDiffuseRadiance.front().radiance,
        cpu.film.indirectDiffuseRadiance.front().radiance,
        "Imported CPU/GPU indirect-diffuse zero plane differs.");
    expectNear(
        gpu.film.indirectSpecularRadiance.front().radiance,
        cpu.film.indirectSpecularRadiance.front().radiance,
        "Imported CPU/GPU indirect-specular zero plane differs.");
    expectCleanValidation(renderer.validationReport(), "Imported CPU/GPU path renderer");
}

int runTest() {
    const std::filesystem::path shaderPath{RAYM0NADE_PATH_TRACE_SHADER};
    VulkanRayQueryOptions vulkanOptions;
    vulkanOptions.requestValidation = true;

    testDeterminismAndPartitionInvariance(shaderPath, vulkanOptions);
    testEnvironmentMissAndEmissiveHit(shaderPath, vulkanOptions);
    testMaterialTexturesCutoutAndTransmission(shaderPath, vulkanOptions);
    testImportedCpuReference(shaderPath, vulkanOptions);
    std::cout << std::fixed << std::setprecision(3)
              << "GPU path-render determinism, partitioning, lighting, Film, and validation tests passed.\n";
    return 0;
}

}  // namespace

int main() {
    try {
        return runTest();
    } catch (const std::exception& error) {
        const std::string_view message{error.what()};
        if (message.rfind(kNoCompatibleDevice, 0U) == 0U) {
            std::cout << "SKIPPED: " << message << '\n';
            return 77;
        }
        std::cerr << "FAILED: " << message << '\n';
        return 1;
    }
}
