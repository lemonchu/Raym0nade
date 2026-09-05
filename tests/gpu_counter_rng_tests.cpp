#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "raym0nade/counter_rng.hpp"
#include "raym0nade/gpu/vulkan_counter_rng_test.hpp"

#ifndef RAYM0NADE_COUNTER_RNG_SHADER
#error "RAYM0NADE_COUNTER_RNG_SHADER must be defined by CMake."
#endif

namespace {

using namespace raym0nade;
using namespace raym0nade::gpu;

constexpr std::string_view kNoCompatibleDevice =
    "No AMD Vulkan device satisfies the Ray Query backend requirements.";

constexpr std::array<CounterRngKey, kVulkanCounterRngKatAddressCount> kKatKeys{{
    {0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U},
    {0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000001U},
    {0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000002U},
    {0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000003U},
    {0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000004U},
    {0x00000001U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U},
    {0x12345678U, 0x9ABCDEF0U, 0x13579BDFU, 0x2468ACE0U, 0x00000000U},
    {0x12345678U, 0x9ABCDEF0U, 0x13579BDFU, 0x2468ACE0U, 0x00000003U},
    {0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU},
    {0x0000002AU, 0x001FA3FFU, 0x0000003FU, 0x00000010U, 0x0000001BU},
    {0xDEADBEEFU, 0x00000011U, 0x0000001DU, 0x00000003U, 0x00000400U},
    {0x00000007U, 0x0001E240U, 0x000F423FU, 0x0000000FU, 0x00000009U},
}};

[[nodiscard]] std::uint32_t floatBits(const float value) noexcept {
    static_assert(sizeof(value) == sizeof(std::uint32_t));
    std::uint32_t result = 0U;
    std::memcpy(&result, &value, sizeof(result));
    return result;
}

void expect(const bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void validateDispatch(
    const std::array<VulkanCounterRngObservation, kVulkanCounterRngKatAddressCount>&
        observations,
    const std::string& description) {
    for (std::size_t index = 0U; index < kKatKeys.size(); ++index) {
        const std::uint32_t expectedWord = counterRandomUint32(kKatKeys[index]);
        const std::uint32_t expectedOpen01Bits =
            floatBits(counterRandomWordToOpen01(expectedWord));
        expect(
            observations[index].word == expectedWord,
            description + " raw word differs from the CPU contract at address " +
                std::to_string(index) + '.');
        expect(
            observations[index].open01Bits == expectedOpen01Bits,
            description + " open-(0,1) float bits differ from the CPU contract at address " +
                std::to_string(index) + '.');
    }
}

int runTest() {
    VulkanRayQueryOptions options;
    options.requestValidation = true;
    const VulkanCounterRngKatResult result = runVulkanCounterRngKat(
        std::filesystem::path{RAYM0NADE_COUNTER_RNG_SHADER}, options);

    expect(!result.deviceName.empty(), "The Vulkan counter RNG KAT has no device name.");
    validateDispatch(result.firstDispatch, "First GPU dispatch");
    validateDispatch(result.repeatedDispatch, "Repeated GPU dispatch");
    for (std::size_t index = 0U; index < result.firstDispatch.size(); ++index) {
        expect(
            result.firstDispatch[index].word == result.repeatedDispatch[index].word &&
                result.firstDispatch[index].open01Bits ==
                    result.repeatedDispatch[index].open01Bits,
            "Repeated GPU dispatch changed address " + std::to_string(index) + '.');
    }

    expect(result.validation.requested, "The GPU KAT did not request Vulkan validation.");
    expect(
        result.validation.errorCount == 0U && result.validation.warningCount == 0U,
        "The GPU KAT reported " + std::to_string(result.validation.errorCount) +
            " Vulkan validation error(s) and " +
            std::to_string(result.validation.warningCount) + " warning(s).");

    std::cout << "GPU counter RNG KAT device: " << result.deviceName << '\n'
              << "Validation: layer="
              << (result.validation.enabled ? "enabled" : "unavailable")
              << ", synchronization="
              << (result.validation.synchronizationValidationEnabled ? "enabled"
                                                                      : "unavailable")
              << " (" << result.validation.errorCount << " error(s), "
              << result.validation.warningCount << " warning(s))\n"
              << "All 12 CPU/GLSL counter RNG addresses matched across two dispatches.\n";
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
