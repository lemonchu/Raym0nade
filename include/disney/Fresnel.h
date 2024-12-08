#ifndef RAYM0NADE_FRESNEL_H
#define RAYM0NADE_FRESNEL_H

#include <glm/glm.hpp>
#include <cmath>
#include <algorithm>

namespace Fresnel {

    glm::vec3 Schlick(const glm::vec3& r0, float radians);

    float Schlick(float r0, float radians);
    float SchlickWeight(float u);
    float SchlickDielectric(float cosThetaI, float relativeIor);
    float Dielectric(float cosThetaI, float ni, float nt);

    float SchlickR0FromRelativeIOR(float eta);

}

#endif //RAYM0NADE_FRESNEL_H