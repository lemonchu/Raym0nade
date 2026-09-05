#pragma once

#include <array>
#include <cstdint>

namespace raym0nade {

// A stable, stateless random key for renderer work. Changing the field-to-counter
// mapping or the algorithm changes the renderer's deterministic sample sequence.
struct CounterRngKey {
    std::uint32_t seed{0U};
    std::uint32_t pixelIndex{0U};
    std::uint32_t sampleIndex{0U};
    std::uint32_t bounceIndex{0U};
    std::uint32_t dimension{0U};
};

using CounterRngBlock = std::array<std::uint32_t, 4>;

namespace detail {

// Exposed for cross-language known-answer tests; renderer code should use the
// address-based APIs below.
[[nodiscard]] CounterRngBlock philox4x32TenRounds(
    const CounterRngBlock& counter,
    const std::array<std::uint32_t, 2>& key) noexcept;

}  // namespace detail

// Returns four adjacent dimensions. dimensionBlock is dimension / 4.
[[nodiscard]] CounterRngBlock counterRandomBlock(
    std::uint32_t seed,
    std::uint32_t pixelIndex,
    std::uint32_t sampleIndex,
    std::uint32_t bounceIndex,
    std::uint32_t dimensionBlock) noexcept;

// Philox4x32-10 maps each key to one deterministic 32-bit word without mutable state.
[[nodiscard]] std::uint32_t counterRandomUint32(const CounterRngKey& key) noexcept;

// Uses the high 23 random bits to select the midpoint of one of 2^23 equal bins.
// The result is always finite and strictly inside the open interval (0, 1).
[[nodiscard]] float counterRandomWordToOpen01(std::uint32_t word) noexcept;

[[nodiscard]] float counterRandomOpen01(const CounterRngKey& key) noexcept;

}  // namespace raym0nade
