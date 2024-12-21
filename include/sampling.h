#ifndef SAMPLING_H
#define SAMPLING_H

#include <random>
#include "geometry.h"
#include "model.h"

class BSDF {
private:
    void sampleGTR2(Generator &gen, vec3 &outDir, vec3 &brdfPdf, float &pdf, int &fails) const;
    void sampleCos(Generator &gen, vec3 &Dir, vec3 &brdfPdf, float &pdf, int &fails) const;
    [[nodiscard]] vec3 getBRDF(vec3 L) const;
    [[nodiscard]] vec3 getBTDF(vec3 L) const;
public:
    vec3 inDir;
    HitInfo surface;

    BSDF(const vec3 &inDir, const vec3 &hit_position);
    BSDF(const vec3 &inDir, const HitInfo &hitInfo);

    [[nodiscard]] vec3 getBSDF(vec3 outDir) const;
    void preciseRefraction(vec3 &outDir, float &F) const;
    void sampleReflection(Generator &gen, vec3 &Dir, vec3 &bsdfPdf, int &fails) const;
};

struct LightSample {
    vec3 bsdfPdf, light;
    float weight;
    LightSample(const vec3 &bsdfPdf, const vec3 &light, float weight);
};

const float eps_lightRadius = 5e-3f;
std::vector<LightSample> sampleDirectLight(const BSDF &bsdf, const Model &model, Generator &gen, int sampleCnt);

#endif // SAMPLING_H