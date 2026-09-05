#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace raym0nade::gpu {

struct VulkanRayQueryObservation {
    bool hit{false};
    std::uint32_t primitiveId{0xffffffffU};
    float distance{0.0F};
    float barycentricU{0.0F};
    float barycentricV{0.0F};
};

struct VulkanRayQuerySelfTestResult {
    std::string deviceName;
    std::string failureReason;
    std::vector<std::string> validationMessages;
    VulkanRayQueryObservation hitRay;
    VulkanRayQueryObservation missRay;
    std::uint32_t validationErrorCount{0};
    std::uint32_t validationWarningCount{0};
    bool validationRequested{false};
    bool validationEnabled{false};
    bool synchronizationValidationEnabled{false};
    bool passed{false};
};

// Runs a headless, deterministic hardware intersection test on a compatible AMD device.
// The SPIR-V file must contain the ray_query_triangle.comp compute shader.
[[nodiscard]] VulkanRayQuerySelfTestResult runVulkanRayQuerySelfTest(
    const std::filesystem::path& spirvPath, bool requestValidation = false);

[[nodiscard]] std::string formatVulkanRayQuerySelfTestResult(
    const VulkanRayQuerySelfTestResult& result);

}  // namespace raym0nade::gpu
