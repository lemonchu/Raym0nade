#include "gen.h"

void RandomNumberGenerator::Init(const std::vector<float>& distribution) {
    prefixSums.resize(distribution.size());
    prefixSums[0] = distribution[0];
    for (size_t i = 1; i < distribution.size(); ++i) {
        prefixSums[i] = prefixSums[i - 1] + distribution[i];
    }
}

int RandomNumberGenerator::operator()(std::mt19937 &gen) const {
    std::uniform_real_distribution<float> dist(0.0f, prefixSums.back());
    float randomValue = dist(gen);
    return std::lower_bound(prefixSums.begin(), prefixSums.end(), randomValue) - prefixSums.begin();
}