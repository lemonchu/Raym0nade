#ifndef GEN_H
#define GEN_H

#include <random>
#include "geometry.h"

class RandomNumberGenerator {
private:
    std::vector<float> prefixSums;
public:
    void Init(const std::vector<float>& distribution);
    int operator()(std::mt19937 &gen) const;
};

class RandomRayGenerator {
private:
    std::vector<vec3> rays;
public:
    RandomRayGenerator(int samples);
    vec3 RandomRayGenerator::operator()(std::mt19937 &gen) const;
};

#endif // GEN_H