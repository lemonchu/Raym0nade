#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <type_traits>

#include "raym0nade/gpu/vulkan_ray_query.hpp"

namespace raym0nade::gpu {

inline constexpr std::size_t kVulkanCounterRngKatAddressCount = 12U;

// Fixed std430 record written by counter_rng_compile_test.comp.
struct alignas(8) VulkanCounterRngObservation {
    std::uint32_t word{0U};
    std::uint32_t open01Bits{0U};
};

struct VulkanCounterRngKatResult {
    std::string deviceName;
    std::array<VulkanCounterRngObservation, kVulkanCounterRngKatAddressCount>
        firstDispatch{};
    std::array<VulkanCounterRngObservation, kVulkanCounterRngKatAddressCount>
        repeatedDispatch{};
    VulkanValidationReport validation;
};

// Dispatches the fixed renderer-address KAT twice on a compatible AMD Vulkan device.
// The SPIR-V file must contain counter_rng_compile_test.comp.
[[nodiscard]] VulkanCounterRngKatResult runVulkanCounterRngKat(
    const std::filesystem::path& spirvPath,
    VulkanRayQueryOptions options = {});

static_assert(sizeof(VulkanCounterRngObservation) == 8U);
static_assert(alignof(VulkanCounterRngObservation) == 8U);
static_assert(offsetof(VulkanCounterRngObservation, word) == 0U);
static_assert(offsetof(VulkanCounterRngObservation, open01Bits) == 4U);
static_assert(std::is_standard_layout_v<VulkanCounterRngObservation>);
static_assert(std::is_trivially_copyable_v<VulkanCounterRngObservation>);

}  // namespace raym0nade::gpu
