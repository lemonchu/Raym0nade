#include "raym0nade/counter_rng.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace raym0nade {
namespace {

constexpr std::uint32_t kPhiloxMultiplier0 = 0xD2511F53U;
constexpr std::uint32_t kPhiloxMultiplier1 = 0xCD9E8D57U;
constexpr std::uint32_t kPhiloxWeyl0 = 0x9E3779B9U;
constexpr std::uint32_t kPhiloxWeyl1 = 0xBB67AE85U;
constexpr int kPhiloxRoundCount = 10;
static_assert(std::numeric_limits<float>::is_iec559);
static_assert(std::numeric_limits<float>::digits == 24);

struct ProductParts {
    std::uint32_t high;
    std::uint32_t low;
};

[[nodiscard]] ProductParts multiplyHighLow(
    const std::uint32_t lhs, const std::uint32_t rhs) noexcept {
    const std::uint64_t product =
        static_cast<std::uint64_t>(lhs) * static_cast<std::uint64_t>(rhs);
    return ProductParts{
        static_cast<std::uint32_t>(product >> 32U),
        static_cast<std::uint32_t>(product),
    };
}

[[nodiscard]] CounterRngBlock philoxRound(
    const CounterRngBlock& counter,
    const std::uint32_t key0,
    const std::uint32_t key1) noexcept {
    const ProductParts first = multiplyHighLow(kPhiloxMultiplier0, counter[0]);
    const ProductParts second = multiplyHighLow(kPhiloxMultiplier1, counter[2]);
    return CounterRngBlock{
        second.high ^ counter[1] ^ key0,
        second.low,
        first.high ^ counter[3] ^ key1,
        first.low,
    };
}

}  // namespace

namespace detail {

CounterRngBlock philox4x32TenRounds(
    const CounterRngBlock& initialCounter,
    const std::array<std::uint32_t, 2>& initialKey) noexcept {
    CounterRngBlock counter = initialCounter;
    std::uint32_t key0 = initialKey[0];
    std::uint32_t key1 = initialKey[1];
    for (int round = 0; round < kPhiloxRoundCount; ++round) {
        counter = philoxRound(counter, key0, key1);
        key0 += kPhiloxWeyl0;
        key1 += kPhiloxWeyl1;
    }
    return counter;
}

}  // namespace detail

CounterRngBlock counterRandomBlock(
    const std::uint32_t seed,
    const std::uint32_t pixelIndex,
    const std::uint32_t sampleIndex,
    const std::uint32_t bounceIndex,
    const std::uint32_t dimensionBlock) noexcept {
    return detail::philox4x32TenRounds(
        CounterRngBlock{pixelIndex, sampleIndex, bounceIndex, dimensionBlock},
        std::array<std::uint32_t, 2>{seed, 0U});
}

std::uint32_t counterRandomUint32(const CounterRngKey& key) noexcept {
    const CounterRngBlock block = counterRandomBlock(
        key.seed, key.pixelIndex, key.sampleIndex, key.bounceIndex, key.dimension >> 2U);
    return block[static_cast<std::size_t>(key.dimension & 3U)];
}

float counterRandomWordToOpen01(const std::uint32_t word) noexcept {
    constexpr float inverseBinCount = 1.0F / 8388608.0F;
    return (static_cast<float>(word >> 9U) + 0.5F) * inverseBinCount;
}

float counterRandomOpen01(const CounterRngKey& key) noexcept {
    return counterRandomWordToOpen01(counterRandomUint32(key));
}

}  // namespace raym0nade
