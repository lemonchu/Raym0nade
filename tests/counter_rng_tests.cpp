#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <numeric>
#include <string_view>

#include "raym0nade/counter_rng.hpp"

namespace {

using namespace raym0nade;

int failureCount = 0;

void expect(const bool condition, const std::string_view message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        ++failureCount;
    }
}

[[nodiscard]] std::uint32_t floatBits(const float value) noexcept {
    static_assert(sizeof(value) == sizeof(std::uint32_t));
    std::uint32_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

struct GoldenVector {
    CounterRngKey key;
    std::uint32_t expectedWord;
    std::uint32_t expectedFloatBits;
};

void testPhiloxCoreGoldenVectors() {
    struct CoreGoldenVector {
        CounterRngBlock counter;
        std::array<std::uint32_t, 2> key;
        CounterRngBlock expected;
    };
    constexpr std::array<CoreGoldenVector, 3> vectors{{
        {{{0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U}},
         {{0x00000000U, 0x00000000U}},
         {{0x6627E8D5U, 0xE169C58DU, 0xBC57AC4CU, 0x9B00DBD8U}}},
        {{{0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU}},
         {{0xFFFFFFFFU, 0xFFFFFFFFU}},
         {{0x408F276DU, 0x41C83B0EU, 0xA20BC7C6U, 0x6D5451FDU}}},
        {{{0x243F6A88U, 0x85A308D3U, 0x13198A2EU, 0x03707344U}},
         {{0xA4093822U, 0x299F31D0U}},
         {{0xD16CFE09U, 0x94FDCCEBU, 0x5001E420U, 0x24126EA1U}}},
    }};
    for (const CoreGoldenVector& vector : vectors) {
        expect(
            detail::philox4x32TenRounds(vector.counter, vector.key) == vector.expected,
            "Philox4x32-10 must match the upstream Random123 known-answer vector.");
    }
}

void testGoldenVectors() {
    constexpr std::array<GoldenVector, 12> vectors{{
        {{0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U}, 0x6627E8D5U, 0x3ECC4FD2U},
        {{0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000001U}, 0xE169C58DU, 0x3F6169C5U},
        {{0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000002U}, 0xBC57AC4CU, 0x3F3C57ADU},
        {{0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000003U}, 0x9B00DBD8U, 0x3F1B00DBU},
        {{0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000004U}, 0x2DCE73E5U, 0x3E3739CCU},
        {{0x00000001U, 0x00000000U, 0x00000000U, 0x00000000U, 0x00000000U}, 0xE3E80670U, 0x3F63E807U},
        {{0x12345678U, 0x9ABCDEF0U, 0x13579BDFU, 0x2468ACE0U, 0x00000000U}, 0xF50A2354U, 0x3F750A23U},
        {{0x12345678U, 0x9ABCDEF0U, 0x13579BDFU, 0x2468ACE0U, 0x00000003U}, 0x6C46BBBEU, 0x3ED88D76U},
        {{0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU, 0xFFFFFFFFU}, 0xB5D79F4CU, 0x3F35D79FU},
        {{0x0000002AU, 0x001FA3FFU, 0x0000003FU, 0x00000010U, 0x0000001BU}, 0xF24F9B31U, 0x3F724F9BU},
        {{0xDEADBEEFU, 0x00000011U, 0x0000001DU, 0x00000003U, 0x00000400U}, 0xB17247A6U, 0x3F317247U},
        {{0x00000007U, 0x0001E240U, 0x000F423FU, 0x0000000FU, 0x00000009U}, 0xDE317E63U, 0x3F5E317FU},
    }};

    for (const GoldenVector& vector : vectors) {
        const std::uint32_t word = counterRandomUint32(vector.key);
        expect(word == vector.expectedWord, "Counter RNG word must match its golden vector.");
        const CounterRngBlock block = counterRandomBlock(
            vector.key.seed,
            vector.key.pixelIndex,
            vector.key.sampleIndex,
            vector.key.bounceIndex,
            vector.key.dimension >> 2U);
        expect(
            block[vector.key.dimension & 3U] == word,
            "The scalar address API must select the matching block lane.");
        expect(
            floatBits(counterRandomOpen01(vector.key)) == vector.expectedFloatBits,
            "Counter RNG float must match its golden vector.");
    }
}

void testOpenIntervalEndpoints() {
    struct ConversionGoldenVector {
        std::uint32_t word;
        std::uint32_t expectedFloatBits;
    };
    constexpr std::array<ConversionGoldenVector, 4> vectors{{
        {0x00000000U, 0x33800000U},
        {0x00000200U, 0x34400000U},
        {0x7FFFFFFFU, 0x3EFFFFFEU},
        {0xFFFFFFFFU, 0x3F7FFFFFU},
    }};
    for (const ConversionGoldenVector& vector : vectors) {
        expect(
            floatBits(counterRandomWordToOpen01(vector.word)) == vector.expectedFloatBits,
            "Open-interval conversion must match its fixed-point golden value.");
    }

    const float minimum = counterRandomWordToOpen01(vectors.front().word);
    const float maximum = counterRandomWordToOpen01(vectors.back().word);
    expect(minimum > 0.0F && minimum < 1.0F, "The lowest converted value must be inside (0, 1).");
    expect(maximum > 0.0F && maximum < 1.0F, "The highest converted value must be inside (0, 1).");

    for (std::uint32_t dimension = 0U; dimension < 250'000U; ++dimension) {
        const float value = counterRandomOpen01(
            CounterRngKey{0xA5A5A5A5U, dimension * 17U, dimension / 11U, dimension / 997U, dimension});
        expect(std::isfinite(value) && value > 0.0F && value < 1.0F,
               "Every sampled value must be finite and strictly inside (0, 1).");
    }
}

void testScheduleOrderIndependence() {
    constexpr std::size_t keyCount = 4096U;
    std::array<CounterRngKey, keyCount> keys{};
    std::array<std::uint32_t, keyCount> reference{};
    for (std::size_t index = 0; index < keyCount; ++index) {
        const auto value = static_cast<std::uint32_t>(index);
        keys[index] = CounterRngKey{
            0xC001D00DU,
            value * 2654435761U,
            value % 257U,
            value % 17U,
            value * 13U,
        };
        reference[index] = counterRandomUint32(keys[index]);
    }

    std::array<std::size_t, keyCount> order{};
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::reverse(order.begin(), order.end());
    std::array<std::uint32_t, keyCount> reverseResults{};
    for (const std::size_t index : order) {
        reverseResults[index] = counterRandomUint32(keys[index]);
    }
    expect(reverseResults == reference, "Reverse scheduling must not change keyed random values.");

    std::array<std::uint32_t, keyCount> permutedResults{};
    for (std::size_t visit = 0; visit < keyCount; ++visit) {
        const std::size_t index = (visit * 4051U + 17U) % keyCount;
        permutedResults[index] = counterRandomUint32(keys[index]);
    }
    expect(permutedResults == reference, "Permuted scheduling must not change keyed random values.");
}

}  // namespace

int main() {
    testPhiloxCoreGoldenVectors();
    testGoldenVectors();
    testOpenIntervalEndpoints();
    testScheduleOrderIndependence();

    if (failureCount == 0) {
        std::cout << "All counter RNG tests passed.\n";
    }
    return failureCount == 0 ? 0 : 1;
}
