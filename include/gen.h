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

#endif // GEN_H