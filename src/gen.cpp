#include <corecrt_math_defines.h>
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

RandomRayGenerator::RandomRayGenerator(int samples) {
    int n = sqrt(samples);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            float u = (i + 0.5f) / n, v = (j + 0.5f) / n * 2.0f * M_PI;
            float d = sqrt( 1.0f - u);
            float z = sqrt(1 - z*z);
            float x = d * cos(v), y = d * sin(v);
            rays.emplace_back(x, y, z);
        }
    }
}

vec3 RandomRayGenerator::operator()(std::mt19937 &gen) const {
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float u = dist(gen), v = dist(gen) * 2.0f * M_PI;
    float d = sqrt( 1.0f - u);
    float z = sqrt(1 - z*z);
    float x = d * cos(v), y = d * sin(v);
    return vec3(x, y, z);
}