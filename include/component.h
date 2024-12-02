#ifndef COMPONENT_H
#define COMPONENT_H

#include <random>
#include "geometry.h"
#include "material.h"
#include "sobel.h"

struct LightFace {
    glm::vec3 position, normal;
    float power;
    LightFace(glm::vec3 position, glm::vec3 normal, float power);
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
    glm::vec3 center, color;
    float power, powerDensity;
    std::vector<LightFace> lightFaces;
    RandomDistribution faceDist;
    LightObject();
};

struct Face {
    vec3 v[3];
    vec2 uv[3];
    Material *material;
    LightObject *lightObject;

    [[nodiscard]] vec3 center() const;
    [[nodiscard]] Box aabb() const;
};

struct HitRecord {
    float t_min, t_max;
    const Face* face;
    HitRecord();
    HitRecord(float t_min, float t_max);
};

#endif // COMPONENT_H