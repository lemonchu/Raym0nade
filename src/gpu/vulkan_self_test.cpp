#include "raym0nade/gpu/vulkan_self_test.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "raym0nade/gpu/vulkan_ray_query.hpp"
#include "raym0nade/scene_data.hpp"

namespace raym0nade::gpu {
namespace {

constexpr float kFloatTolerance = 1.0e-4F;

PackedSceneData makeSelfTestScene() {
    PackedSceneData scene;
    scene.vertices = {
        PackedVertex{{-1.0F, -1.0F, 0.0F, 0.0F}, {0.0F, 1.0F, 0.0F, 0.0F}},
        PackedVertex{{1.0F, -1.0F, 0.0F, 0.0F}, {0.0F, 1.0F, 1.0F, 0.0F}},
        PackedVertex{{0.0F, 1.0F, 0.0F, 0.0F}, {0.0F, 1.0F, 0.5F, 1.0F}},
        PackedVertex{{3.0F, -1.0F, 0.0F, 0.0F}, {0.0F, 1.0F, 0.0F, 0.0F}},
        PackedVertex{{5.0F, -1.0F, 0.0F, 0.0F}, {0.0F, 1.0F, 1.0F, 0.0F}},
        PackedVertex{{4.0F, 1.0F, 0.0F, 0.0F}, {0.0F, 1.0F, 0.5F, 1.0F}},
    };
    scene.triangleIndices = {0U, 1U, 2U, 3U, 4U, 5U};
    scene.triangleMaterialIds = {0U, 0U};
    scene.materials.emplace_back();
    return scene;
}

std::vector<VulkanRayQueryRay> makeSelfTestRays() {
    return {
        VulkanRayQueryRay{{0.0F, 0.0F, 1.0F, 0.001F}, {0.0F, 0.0F, -1.0F, 100.0F}},
        VulkanRayQueryRay{{2.0F, 2.0F, 1.0F, 0.001F}, {0.0F, 0.0F, -1.0F, 100.0F}},
    };
}

VulkanRayQueryObservation toObservation(const VulkanRayQueryHit& hit) noexcept {
    VulkanRayQueryObservation result;
    result.hit = hit.hasHit();
    result.primitiveId = hit.primitiveId;
    result.distance = hit.distance;
    result.barycentricU = hit.barycentricU;
    result.barycentricV = hit.barycentricV;
    return result;
}

bool nearlyEqual(float left, float right) noexcept {
    return std::isfinite(left) && std::abs(left - right) <= kFloatTolerance;
}

std::string validateObservations(const VulkanRayQuerySelfTestResult& result) {
    if (!result.hitRay.hit) {
        return "The deterministic hit ray missed the triangle.";
    }
    if (result.hitRay.primitiveId != 0U) {
        return "The hit ray returned an unexpected primitive ID.";
    }
    if (!nearlyEqual(result.hitRay.distance, 1.0F)) {
        return "The hit ray returned an unexpected intersection distance.";
    }
    if (!nearlyEqual(result.hitRay.barycentricU, 0.25F) ||
        !nearlyEqual(result.hitRay.barycentricV, 0.5F)) {
        return "The hit ray returned unexpected barycentric coordinates.";
    }
    if (result.missRay.hit) {
        return "The deterministic miss ray unexpectedly hit the triangle.";
    }
    if (result.missRay.primitiveId != kInvalidRayQueryPrimitiveId) {
        return "The miss ray did not preserve the invalid primitive sentinel.";
    }
    return {};
}

void copyValidationReport(
    const VulkanValidationReport& source, VulkanRayQuerySelfTestResult& destination) {
    destination.validationEnabled = source.enabled;
    destination.synchronizationValidationEnabled =
        source.synchronizationValidationEnabled;
    destination.validationErrorCount = source.errorCount;
    destination.validationWarningCount = source.warningCount;
    destination.validationMessages = source.messages;
}

}  // namespace

VulkanRayQuerySelfTestResult runVulkanRayQuerySelfTest(
    const std::filesystem::path& spirvPath, bool requestValidation) {
    VulkanRayQuerySelfTestResult result;
    result.validationRequested = requestValidation;
    try {
        VulkanRayQueryOptions options;
        options.requestValidation = requestValidation;
        VulkanRayQueryIntersector intersector{makeSelfTestScene(), spirvPath, options};
        result.deviceName = intersector.deviceName();
        const std::vector<VulkanRayQueryRay> rays = makeSelfTestRays();
        const std::vector<VulkanRayQueryRay> firstBatch{
            VulkanRayQueryRay{{4.0F, 0.0F, 1.0F, 0.001F},
                              {0.0F, 0.0F, -1.0F, 100.0F}}};
        const VulkanRayQueryBatch firstResult = intersector.intersect(firstBatch);
        if (firstResult.hits.size() != 1U || !firstResult.hits.front().hasHit() ||
            firstResult.hits.front().primitiveId != 1U) {
            throw std::runtime_error("The persistent Vulkan intersector failed its first batch.");
        }
        const VulkanRayQueryBatch batch = intersector.intersect(rays);
        if (batch.hits.size() != 2U) {
            throw std::runtime_error("The Vulkan self-test returned an unexpected hit count.");
        }
        result.hitRay = toObservation(batch.hits[0]);
        result.missRay = toObservation(batch.hits[1]);
        copyValidationReport(intersector.validationReport(), result);
        result.failureReason = validateObservations(result);
        if (result.failureReason.empty() && result.validationErrorCount != 0U) {
            result.failureReason = "Vulkan validation reported " +
                                   std::to_string(result.validationErrorCount) + " error(s).";
        }
        result.passed = result.failureReason.empty();
    } catch (const std::exception& error) {
        result.failureReason = error.what();
        result.passed = false;
    }
    return result;
}

std::string formatVulkanRayQuerySelfTestResult(const VulkanRayQuerySelfTestResult& result) {
    std::ostringstream output;
    output << "Vulkan Ray Query self-test: " << (result.passed ? "PASS" : "FAIL") << '\n';
    if (!result.deviceName.empty()) {
        output << "Device: " << result.deviceName << '\n';
    }
    if (result.validationRequested) {
        output << "Validation: " << (result.validationEnabled ? "enabled" : "unavailable")
               << ", synchronization="
               << (result.synchronizationValidationEnabled ? "enabled" : "unavailable")
               << " (" << result.validationErrorCount << " error(s), "
               << result.validationWarningCount << " warning(s))\n";
    } else {
        output << "Validation: not requested\n";
    }
    output << std::setprecision(7);
    output << "Hit ray: hit=" << (result.hitRay.hit ? "true" : "false")
           << ", primitive=" << result.hitRay.primitiveId
           << ", t=" << result.hitRay.distance
           << ", barycentrics=(" << result.hitRay.barycentricU << ", "
           << result.hitRay.barycentricV << ")\n";
    output << "Miss ray: hit=" << (result.missRay.hit ? "true" : "false")
           << ", primitive=" << result.missRay.primitiveId << '\n';
    if (!result.failureReason.empty()) {
        output << "Failure: " << result.failureReason << '\n';
    }
    for (const std::string& message : result.validationMessages) {
        output << "Validation message: " << message << '\n';
    }
    return output.str();
}

}  // namespace raym0nade::gpu
