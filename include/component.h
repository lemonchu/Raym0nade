#ifndef COMPONENT_H
#define COMPONENT_H

#include <random>
#include "geometry.h"
#include "material.h"

struct VertexData {
    vec2 uv;
    vec3 normal;

    VertexData();
    VertexData(const vec2 &uv, const vec3 &normal);
};

struct Face {
    vec3 v[3];
    VertexData *data[3];
    const Material *material;

    [[nodiscard]] vec3 center() const;
    [[nodiscard]] Box aabb() const;
};

struct HitRecord {
    float t_min, t_max;
    const Face* face;
    HitRecord();
    HitRecord(float t_min, float t_max);
};

struct Generator {
    std::mt19937 mt;
    std::uniform_real_distribution<float> U;
    explicit Generator(unsigned int seed);
    float operator()();
};

class RandomDistribution {
private:
    std::vector<float> prefixSums;
public:
    void Init(const std::vector<float>& distribution);
    int operator()(Generator &gen) const;
};

class LightObject {
public:
    vec3 center, color;
    float power, powerDensity;
    std::vector<Face> faces;
    RandomDistribution faceDist;
    LightObject();
};

#endif // COMPONENT_H