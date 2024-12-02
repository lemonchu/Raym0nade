#ifndef SOBEL_H
#define SOBEL_H

#include <cstdint>
#include <random>
#include <functional>

class SobelGenerator {
    static const int bufferSize = 128;
    unsigned int cur, buffer[bufferSize];
    unsigned int get(unsigned int i);
public:
    static const int bitDepth = 32;
    unsigned int C[bitDepth];
    SobelGenerator();
    void gen();
    float operator()();
};

class SobelGroup {
public:
    static const int groupSize = 16;
private:
    SobelGenerator gen[groupSize];
    std::mt19937 mt;
    std::uniform_real_distribution<float> U;
public:
    SobelGroup();
    float operator ()();
    float operator ()(int dimension);
};

struct Generator {
    SobelGroup *sobel;
    int dimension;
    float operator ()();
};

#endif // SOBEL_H