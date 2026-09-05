#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace raym0nade::gpu {

struct VulkanVersion {
    std::uint32_t major{0};
    std::uint32_t minor{0};
    std::uint32_t patch{0};
};

struct VulkanDeviceCapabilities {
    std::string deviceName;
    std::string driverName;
    std::string driverInfo;
    std::vector<std::string> missingRequirements;
    VulkanVersion apiVersion;
    std::uint64_t deviceLocalMemoryBytes{0};
    std::uint64_t maximumAccelerationStructurePrimitiveCount{0};
    std::uint32_t vendorId{0};
    std::uint32_t deviceId{0};
    std::uint32_t subgroupSize{0};
    std::uint32_t computeQueueFamily{0};
    bool hasComputeQueue{false};
    bool integrated{false};
    bool bufferDeviceAddress{false};
    bool accelerationStructure{false};
    bool rayQuery{false};
    bool rayTracingPipeline{false};

    [[nodiscard]] bool isAmd() const noexcept;
    [[nodiscard]] bool supportsRayQueryBackend() const noexcept;
};

[[nodiscard]] VulkanVersion vulkanLoaderVersion();
[[nodiscard]] std::vector<VulkanDeviceCapabilities> enumerateVulkanDevices();

}  // namespace raym0nade::gpu
