#ifndef SAMPLING_H
#define SAMPLING_H

#include <random>
#include "geometry.h"
#include "model.h"

class BRDF {
public:
    vec3 inDir, shapeNormal, surfaceNormal, tangent, bitangent;
    float max, EP_accept, roughness;
    BRDF(const vec3 &inDir, const vec3 &shapeNormal, const vec3 &surfaceNormal);
    vec3 sample(std::mt19937 &gen, float &P_success) const;
    [[nodiscard]] float pdf(const vec3 &outDir) const;
    [[nodiscard]] float P_accept(const vec3 &outDir) const;
};

const int maxTrys = 16;
const float eps_lightRadius = 1e-5f;

vec3 sampleDirectLight(const vec3 &pos, const BRDF &brdf, const Model &model, std::mt19937 &gen);

#endif // SAMPLING_H