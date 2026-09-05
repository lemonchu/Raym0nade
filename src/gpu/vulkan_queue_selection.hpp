#pragma once

#include <cstddef>
#include <cstdint>

namespace raym0nade::gpu::detail {

struct ComputeQueueFamilyInfo {
    std::uint32_t family{0U};
    std::uint32_t queueCount{0U};
    bool supportsCompute{false};
    bool supportsGraphics{false};
};

struct ComputeQueueFamilySelection {
    bool available{false};
    std::uint32_t family{0U};
    std::uint32_t queueCount{0U};
};

[[nodiscard]] constexpr ComputeQueueFamilySelection selectComputeQueueFamily(
    const ComputeQueueFamilyInfo* families,
    const std::size_t familyCount,
    const std::uint32_t requiredQueueCount) noexcept {
    if (requiredQueueCount == 0U) {
        return {};
    }
    for (std::size_t index = 0U; index < familyCount; ++index) {
        const ComputeQueueFamilyInfo& family = families[index];
        if (family.supportsCompute && !family.supportsGraphics &&
            family.queueCount >= requiredQueueCount) {
            return {true, family.family, family.queueCount};
        }
    }
    for (std::size_t index = 0U; index < familyCount; ++index) {
        const ComputeQueueFamilyInfo& family = families[index];
        if (family.supportsCompute && family.queueCount >= requiredQueueCount) {
            return {true, family.family, family.queueCount};
        }
    }
    return {};
}

[[nodiscard]] constexpr ComputeQueueFamilySelection selectLargestComputeQueueFamily(
    const ComputeQueueFamilyInfo* families,
    const std::size_t familyCount) noexcept {
    std::uint32_t maximumQueueCount = 0U;
    for (std::size_t index = 0U; index < familyCount; ++index) {
        const ComputeQueueFamilyInfo& family = families[index];
        if (family.supportsCompute && family.queueCount > maximumQueueCount) {
            maximumQueueCount = family.queueCount;
        }
    }
    return selectComputeQueueFamily(families, familyCount, maximumQueueCount);
}

}  // namespace raym0nade::gpu::detail
