#ifndef SAMPLING_H
#define SAMPLING_H

#include <random>
#include "geometry.h"
#include "model.h"

class BRDF {
public:
    vec3 inDir, tangent, bitangent;
    glm::mat3 worldToTangent;
    HitInfo surface;

    BRDF(const vec3 &inDir);

    void genTangentSpace();
    vec3 getBRDF(vec3 outDir) const;
    vec3 getBRDFWithCos(vec3 outDir) const;
    void sample(Generator &gen, vec3 &Dir, vec3 &brdfPdf, int &trys) const;
};

const float eps_lightRadius = 1e-5f;

vec3 sampleDirectLight(const vec3 &pos, const BRDF &brdf, const Model &model,
                       Generator &gen, vec3 &brdfPdf, int &trys);

#endif // SAMPLING_H