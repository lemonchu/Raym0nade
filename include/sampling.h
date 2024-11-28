#ifndef SAMPLING_H
#define SAMPLING_H

#include <random>
#include "geometry.h"
#include "model.h"

class BRDF {
    vec3 inDir, shapeNormal, surfaceNormal, tangent, bitangent;
    float max, roughness;
public:
    BRDF(const vec3 &inDir, const vec3 &shapeNormal, const vec3 &surfaceNormal);
    vec3 sample(std::mt19937 &gen) const;
    float pdf(const vec3 &outDir) const;
    float P_accept(const vec3 &outDir) const;
};

const int maxTrys = 16;
const float eps_lightRadius = 1e-3f;

void sampleLightObject(const vec3 &pos, const BRDF &brdf, const Model &model, std::mt19937 &gen, float &prob, int &lightIndex);

void sampleLightFace(const vec3 &pos, const LightObject &lightObject, std::mt19937 &gen, int &failCount, int &faceIndex);

void sampleDirectLight(const vec3 &pos, const BRDF &brdf, const Model &model, std::mt19937 &gen, int &failCount, vec3 &light);

const int probeSize = 32;

class Probe {
public:
    int sampleCount;
    vec3 normal, tangent, bitangent;
    float powerSum[probeSize*probeSize], Epower;
    float power(int index) const;
    vec3 sample(std::mt19937 &gen, const BRDF &brdf, float &prob) const;
    Probe();
};

#endif // SAMPLING_H