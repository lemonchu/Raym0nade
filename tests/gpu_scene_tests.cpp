#include <array>
#include <chrono>
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

#include "raym0nade/gpu/vulkan_ray_query.hpp"
#include "raym0nade/model.hpp"
#include "raym0nade/scene_data.hpp"

#ifndef RAYM0NADE_TEST_SOURCE_DIR
#error "RAYM0NADE_TEST_SOURCE_DIR must be defined by CMake."
#endif

#ifndef RAYM0NADE_RAY_QUERY_SCENE_SHADER
#error "RAYM0NADE_RAY_QUERY_SCENE_SHADER must be defined by CMake."
#endif

namespace {

using namespace raym0nade;
using namespace raym0nade::gpu;

constexpr float kComparisonTolerance = 1.0e-4F;
constexpr float kRayMaximum = 10.0F;
constexpr std::string_view kNoCompatibleDevice =
    "No AMD Vulkan device satisfies the Ray Query backend requirements.";

[[noreturn]] void fail(std::size_t rayIndex, const std::string& message) {
    throw std::runtime_error(
        "Ray " + std::to_string(rayIndex) + ": " + message);
}

bool nearlyEqual(float left, float right) noexcept {
    return std::isfinite(left) && std::isfinite(right) &&
           std::abs(left - right) <= kComparisonTolerance;
}

vec3 packedPosition(const PackedVertex& vertex) noexcept {
    return vec3{
        vertex.positionAndNormalX[0],
        vertex.positionAndNormalX[1],
        vertex.positionAndNormalX[2],
    };
}

VulkanRayQueryRay makeGpuRay(const Ray& ray) noexcept {
    return VulkanRayQueryRay{
        {ray.origin.x, ray.origin.y, ray.origin.z, kRayEpsilon},
        {ray.direction.x, ray.direction.y, ray.direction.z, kRayMaximum},
    };
}

void validateGpuHitFields(const VulkanRayQueryHit& hit, std::size_t rayIndex) {
    if (!std::isfinite(hit.distance) || !std::isfinite(hit.barycentricU) ||
        !std::isfinite(hit.barycentricV)) {
        fail(rayIndex, "the GPU returned non-finite intersection fields.");
    }
    if (!hit.hasHit()) {
        return;
    }
    if (hit.distance + kComparisonTolerance < kRayEpsilon ||
        hit.distance > kRayMaximum + kComparisonTolerance) {
        fail(rayIndex, "the GPU intersection distance is outside the requested interval.");
    }
    if (hit.barycentricU < -kComparisonTolerance ||
        hit.barycentricV < -kComparisonTolerance ||
        hit.barycentricU > 1.0F + kComparisonTolerance ||
        hit.barycentricV > 1.0F + kComparisonTolerance ||
        hit.barycentricU + hit.barycentricV > 1.0F + kComparisonTolerance) {
        fail(rayIndex, "the GPU returned out-of-range barycentric coordinates.");
    }
}

void validatePackedPrimitive(
    const PackedSceneData& scene, const HitRecord& cpuHit, std::size_t rayIndex) {
    if (cpuHit.face == nullptr) {
        fail(rayIndex, "a hit has no CPU face pointer.");
    }
    if (cpuHit.primitiveIndex >= scene.triangleCount() ||
        cpuHit.primitiveIndex > std::numeric_limits<std::uint32_t>::max()) {
        fail(rayIndex, "the CPU primitive index is outside the packed scene range.");
    }

    const std::size_t indexOffset = cpuHit.primitiveIndex * 3U;
    for (std::size_t corner = 0U; corner < 3U; ++corner) {
        const std::uint32_t vertexIndex = scene.triangleIndices[indexOffset + corner];
        if (vertexIndex >= scene.vertices.size()) {
            fail(rayIndex, "the packed primitive contains an invalid vertex index.");
        }
        const vec3 position = packedPosition(scene.vertices[vertexIndex]);
        const vec3& cpuPosition = cpuHit.face->vertices[corner];
        if (!nearlyEqual(position.x, cpuPosition.x) ||
            !nearlyEqual(position.y, cpuPosition.y) ||
            !nearlyEqual(position.z, cpuPosition.z)) {
            fail(rayIndex, "the CPU face and packed primitive positions differ.");
        }
    }
}

void validateBatch(
    const std::array<Ray, 5>& rays,
    const std::array<HitRecord, 5>& cpuHits,
    const PackedSceneData& scene,
    const VulkanRayQueryBatch& batch,
    std::size_t batchIndex) {
    if (batch.hits.size() != rays.size()) {
        throw std::runtime_error(
            "GPU batch " + std::to_string(batchIndex) +
            " returned an unexpected result count.");
    }

    for (std::size_t rayIndex = 0U; rayIndex < rays.size(); ++rayIndex) {
        const HitRecord& cpuHit = cpuHits[rayIndex];
        const VulkanRayQueryHit& gpuHit = batch.hits[rayIndex];
        const bool cpuHasHit = cpuHit.face != nullptr;
        validateGpuHitFields(gpuHit, rayIndex);

        if (!cpuHasHit) {
            if (cpuHit.primitiveIndex != std::numeric_limits<std::size_t>::max()) {
                fail(rayIndex, "the CPU miss did not preserve its primitive sentinel.");
            }
            if (gpuHit.hasHit() ||
                gpuHit.primitiveId != kInvalidRayQueryPrimitiveId) {
                fail(rayIndex, "the GPU miss did not preserve both miss sentinels.");
            }
            continue;
        }

        validatePackedPrimitive(scene, cpuHit, rayIndex);
        if (!(cpuHit.tMaximum >= kRayEpsilon && cpuHit.tMaximum <= kRayMaximum)) {
            fail(rayIndex, "the CPU hit lies outside the requested ray interval.");
        }
        if (!gpuHit.hasHit()) {
            fail(rayIndex, "the GPU missed a CPU reference hit.");
        }

        const std::uint32_t expectedPrimitive =
            static_cast<std::uint32_t>(cpuHit.primitiveIndex);
        if (gpuHit.primitiveId != expectedPrimitive) {
            fail(rayIndex, "the CPU and GPU primitive indices differ.");
        }
        if (!nearlyEqual(gpuHit.distance, cpuHit.tMaximum)) {
            fail(rayIndex, "the CPU and GPU intersection distances differ.");
        }

        const vec3 point = rays[rayIndex].origin +
                           cpuHit.tMaximum * rays[rayIndex].direction;
        const vec3 cpuBarycentrics = barycentric(
            cpuHit.face->vertices[0],
            cpuHit.face->vertices[1],
            cpuHit.face->vertices[2],
            point);
        if (!isFinite(cpuBarycentrics) ||
            !nearlyEqual(gpuHit.barycentricU, cpuBarycentrics.y) ||
            !nearlyEqual(gpuHit.barycentricV, cpuBarycentrics.z)) {
            fail(rayIndex, "the CPU and GPU barycentric coordinates differ.");
        }
    }
}

void validateGenericIntersectorRejectsCutout(
    const PackedSceneData& sourceScene,
    const std::filesystem::path& shaderPath) {
    PackedSceneData cutoutScene = sourceScene;
    const std::uint32_t materialId = cutoutScene.triangleMaterialIds.front();
    cutoutScene.textures.push_back(PackedTexture{0U, 1U, 1U, 1U});
    cutoutScene.textureMipLevels.push_back(
        PackedTextureMip{0U, 1U, 1U, 1U});
    cutoutScene.textureTexelsRgba8.push_back(0x00ffffffU);
    PackedMaterial& material = cutoutScene.materials[materialId];
    material.textureIds[0] = 0U;
    material.flagsAndReserved[0] |=
        kPackedMaterialHasDiffuseTexture | kPackedMaterialCutout;
    cutoutScene.validate();

    bool rejected = false;
    try {
        VulkanRayQueryIntersector unsupported{
            cutoutScene, shaderPath, VulkanRayQueryOptions{}};
        (void)unsupported;
    } catch (const std::invalid_argument& error) {
        rejected =
            std::string_view{error.what()}.find("alpha-cutout") !=
            std::string_view::npos;
    }
    if (!rejected) {
        throw std::runtime_error(
            "The generic Vulkan intersector must reject alpha-cutout scenes.");
    }
}

int runTest() {
    const std::filesystem::path sourceDirectory{RAYM0NADE_TEST_SOURCE_DIR};
    Model model{sourceDirectory / "tests" / "data", "triangle.obj", "null"};
    const auto packBegin = std::chrono::steady_clock::now();
    PackedSceneData scene = model.packScene();
    const auto packEnd = std::chrono::steady_clock::now();
    if (scene.areaLights.empty() || scene.areaLightTriangles.empty()) {
        throw std::runtime_error(
            "The GPU scene fixture must exercise non-empty area-light uploads.");
    }
    scene.environment = PackedEnvironment{
        2U,
        1U,
        kPackedEnvironmentHasImportance,
        0U,
    };
    scene.environmentRows = {
        PackedEnvironmentRow{1.0F, 1.0F, 6.28318548F, 0.0F},
    };
    scene.environmentTexels = {
        PackedEnvironmentTexel{{1.0F, 1.0F, 1.0F, 0.5F}},
        PackedEnvironmentTexel{{1.0F, 1.0F, 1.0F, 1.0F}},
    };
    scene.validate();

    const std::array<Ray, 5> rays{
        Ray{vec3{0.0F, 0.0F, 0.0F}, vec3{0.0F, 0.0F, 1.0F}},
        Ray{vec3{-0.5F, -0.5F, 0.0F}, vec3{0.0F, 0.0F, 1.0F}},
        Ray{vec3{1.1F, 1.1F, 0.0F}, vec3{0.0F, 0.0F, 1.0F}},
        Ray{vec3{0.0F, 0.0F, 3.0F}, vec3{0.0F, 0.0F, -1.0F}},
        Ray{vec3{2.0F, 0.0F, 0.0F}, vec3{0.0F, 0.0F, 1.0F}},
    };
    constexpr std::array<bool, 5> kExpectedHits{true, true, true, true, false};

    std::array<HitRecord, 5> cpuHits{};
    std::vector<VulkanRayQueryRay> gpuRays;
    gpuRays.reserve(rays.size());
    for (std::size_t rayIndex = 0U; rayIndex < rays.size(); ++rayIndex) {
        cpuHits[rayIndex] = model.intersect(rays[rayIndex]);
        if ((cpuHits[rayIndex].face != nullptr) != kExpectedHits[rayIndex]) {
            fail(rayIndex, "the imported CPU reference scene produced an unexpected result.");
        }
        gpuRays.push_back(makeGpuRay(rays[rayIndex]));
    }

    VulkanRayQueryOptions options;
    options.requestValidation = true;
    const std::filesystem::path shaderPath{
        RAYM0NADE_RAY_QUERY_SCENE_SHADER};
    validateGenericIntersectorRejectsCutout(scene, shaderPath);
    const auto setupBegin = std::chrono::steady_clock::now();
    VulkanRayQueryIntersector intersector{
        scene, shaderPath, options};
    const auto setupEnd = std::chrono::steady_clock::now();

    const VulkanRayQueryBatch firstBatch = intersector.intersect(gpuRays);
    const VulkanRayQueryBatch secondBatch = intersector.intersect(gpuRays);
    validateBatch(rays, cpuHits, scene, firstBatch, 1U);
    validateBatch(rays, cpuHits, scene, secondBatch, 2U);
    const VulkanValidationReport validation = intersector.validationReport();
    if (validation.errorCount != 0U) {
        throw std::runtime_error(
            "Vulkan validation reported " + std::to_string(validation.errorCount) +
            " error(s).");
    }

    const double packedSceneMilliseconds =
        std::chrono::duration<double, std::milli>{packEnd - packBegin}.count();
    const double setupHostMilliseconds =
        std::chrono::duration<double, std::milli>{setupEnd - setupBegin}.count();
    const VulkanRayQuerySetupTimings& setup = intersector.setupTimings();
    std::cout << std::fixed << std::setprecision(3)
              << "GPU scene test device: " << intersector.deviceName() << '\n'
              << "Packed-scene conversion: " << packedSceneMilliseconds << " ms\n"
              << "Setup host total: " << setupHostMilliseconds << " ms"
              << " (upload " << setup.uploadMilliseconds << " ms, AS build "
              << setup.accelerationStructureBuildMilliseconds << " ms)\n"
              << "Batch 1 host dispatch/readback: "
              << firstBatch.dispatchAndReadbackMilliseconds << " ms\n"
              << "Batch 2 host dispatch/readback: "
              << secondBatch.dispatchAndReadbackMilliseconds << " ms\n"
              << "Validation: requested=" << (validation.requested ? "yes" : "no")
              << ", layer=" << (validation.enabled ? "enabled" : "unavailable")
              << ", synchronization="
              << (validation.synchronizationValidationEnabled ? "enabled" : "unavailable")
              << " (" << validation.errorCount << " error(s), "
              << validation.warningCount << " warning(s))\n"
              << "CPU/GPU packed-scene intersection test passed.\n";
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
