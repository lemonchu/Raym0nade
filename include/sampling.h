#ifndef SAMPLING_H
#define SAMPLING_H

#include <random>
#include "geometry.h"
#include "model.h"

class BRDF {
public:
    vec3 inDir, shapeNormal, surfaceNormal, tangent, bitangent;
    glm::mat3 worldToTangent;
    float max, EP_accept;
    vec3 baseColor;
    vec3 transmittanceColor;

    float sheen, sheenTint;
    float clearcoat, clearcoatGloss;
    float metallic;
    float specTrans, diffTrans;
    float flatness;
    float anisotropic;
    float relativeIOR;
    float specularTint;
    float roughness;
    float scatterDistance;
    float ior;

    BRDF(const vec3 &inDir, const vec3 &shapeNormal, const vec3 &surfaceNormal);

    vec3 sample(Generator &gen, float &P_success) const;
    [[nodiscard]] float pdf(const vec3 &outDir) const;
    [[nodiscard]] float P_accept(const vec3 &outDir) const;
};

const int MaxTrys = 16;
const float eps_lightRadius = 1e-5f;

vec3 sampleDirectLight(const vec3 &pos, const BRDF &brdf, const Model &model, Generator &gen);

#endif // SAMPLING_H