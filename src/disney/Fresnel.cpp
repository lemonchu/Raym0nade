#include "disney/Fresnel.h"

namespace Fresnel {
    glm::vec3 Schlick(const glm::vec3& r0, float radians) {
        float exponential = std::pow(1.0f - radians, 5.0f);
        return r0 + (glm::vec3(1.0f) - r0) * exponential;
    }

    float Schlick(float r0, float radians) {
        return glm::mix(1.0f, Fresnel::SchlickWeight(radians), r0);
    }

    float SchlickWeight(float u) {
        float m = std::clamp(1.0f - u, 0.0f, 1.0f);
        float m2 = m * m;
        return m * m2 * m2;
    }

    float SchlickDielectic(float cosThetaI, float relativeIor) {
        float r0 = SchlickR0FromRelativeIOR(relativeIor);
        return r0 + (1.0f - r0) * SchlickWeight(cosThetaI);
    }

    float Dielectric(float cosThetaI, float ni, float nt) {
        cosThetaI = std::clamp(cosThetaI, -1.0f, 1.0f);

        // Swap index of refraction if this is coming from inside the surface
        if (cosThetaI < 0.0f) {
            std::swap(ni, nt);
            cosThetaI = -cosThetaI;
        }

        float sinThetaI = std::sqrt(std::fmax(0.0f, 1.0f - cosThetaI * cosThetaI));
        float sinThetaT = ni / nt * sinThetaI;

        // Check for total internal reflection
        if (sinThetaT >= 1) {
            return 1;
        }

        float cosThetaT = std::sqrt(std::fmax(0.0f, 1.0f - sinThetaT * sinThetaT));

        float rParallel = ((nt * cosThetaI) - (ni * cosThetaT)) / ((nt * cosThetaI) + (ni * cosThetaT));
        float rPerpendicuar = ((ni * cosThetaI) - (nt * cosThetaT)) / ((ni * cosThetaI) + (nt * cosThetaT));
        return (rParallel * rParallel + rPerpendicuar * rPerpendicuar) / 2;
    }

    float SchlickR0FromRelativeIOR(float eta) {
        return glm::pow(eta - 1.0f, 2.0f) / glm::pow(eta + 1.0f, 2.0f);
    }
}