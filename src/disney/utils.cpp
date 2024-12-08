#include "disney/utils.h"

namespace Bsdf {
    float SeparableSmithGGXG1(const glm::vec3 &w, const glm::vec3 &wm, float ax, float ay) {
        float absTanTheta = std::abs(std::tan(w.y));
        if (std::isinf(absTanTheta)) {
            return 0.0f;
        }

        float a = std::sqrt(std::cos(w.y) * std::cos(w.y) * ax * ax + std::sin(w.y) * std::sin(w.y) * ay * ay);
        float a2Tan2Theta = (a * absTanTheta) * (a * absTanTheta);

        float lambda = 0.5f * (-1.0f + std::sqrt(1.0f + a2Tan2Theta));
        return 1.0f / (1.0f + lambda);
    }

    float SeparableSmithGGXG1(const glm::vec3 &w, float a) {
        float a2 = a * a;
        float absDotNV = std::abs(std::cos(w.y));

        return 2.0f / (1.0f + std::sqrt(a2 + (1 - a2) * absDotNV * absDotNV));
    }

    float HeightCorrelatedSmithGGXG2(const glm::vec3 &wo, const glm::vec3 &wi, float a) {
        float absDotNV = std::abs(std::cos(wo.y));
        float absDotNL = std::abs(std::cos(wi.y));
        float a2 = a * a;

        float denomA = absDotNV * std::sqrt(a2 + (1.0f - a2) * absDotNL * absDotNL);
        float denomB = absDotNL * std::sqrt(a2 + (1.0f - a2) * absDotNV * absDotNV);

        return 2.0f * absDotNL * absDotNV / (denomA + denomB);
    }

    float GgxIsotropicD(const glm::vec3 &wm, float a) {
        float a2 = a * a;
        float dotNH2 = std::cos(wm.y) * std::cos(wm.y);

        float sqrtdenom = dotNH2 * (a2 - 1) + 1;
        return a2 / (M_PI * sqrtdenom * sqrtdenom);
    }

    float GgxAnisotropicD(const glm::vec3 &wm, float ax, float ay) {
        float dotHX2 = wm.x * wm.x;
        float dotHY2 = wm.z * wm.z;
        float cos2Theta = std::cos(wm.y) * std::cos(wm.y);
        float ax2 = ax * ax;
        float ay2 = ay * ay;

        return 1.0f / (M_PI * ax * ay * ((dotHX2 / ax2) + (dotHY2 / ay2) + cos2Theta) *
                       ((dotHX2 / ax2) + (dotHY2 / ay2) + cos2Theta));
    }

    glm::vec3 SampleGgxVndf(glm::vec3 wo, float roughness, float u1, float u2) {
        return SampleGgxVndfAnisotropic(wo, roughness, roughness, u1, u2);
    }

    float GgxVndfPdf(const glm::vec3 &wo, const glm::vec3 &wm, const glm::vec3 &wi, float a) {
        float absDotNL = std::abs(std::cos(wi.y));
        float absDotLH = std::abs(glm::dot(wm, wi));

        float G1 = Bsdf::SeparableSmithGGXG1(wo, a);
        float D = Bsdf::GgxIsotropicD(wm, a);

        return G1 * absDotLH * D / absDotNL;
    }

    glm::vec3 SampleGgxVndfAnisotropic(const glm::vec3 &wo, float ax, float ay, float u1, float u2) {
        glm::vec3 v = glm::normalize(glm::vec3(wo.x * ax, wo.y, wo.z * ay));

        glm::vec3 t1 = (v.y < 0.9999f) ? glm::normalize(glm::cross(v, glm::vec3(0.0f, 1.0f, 0.0f))) : glm::vec3(1.0f,
                                                                                                                0.0f,
                                                                                                                0.0f);
        glm::vec3 t2 = glm::cross(t1, v);

        float a = 1.0f / (1.0f + v.y);
        float r = std::sqrt(u1);
        float phi = (u2 < a) ? (u2 / a) * M_PI : M_PI + (u2 - a) / (1.0f - a) * M_PI;
        float p1 = r * std::cos(phi);
        float p2 = r * std::sin(phi) * ((u2 < a) ? 1.0f : v.y);

        glm::vec3 n = p1 * t1 + p2 * t2 + std::sqrt(std::max(0.0f, 1.0f - p1 * p1 - p2 * p2)) * v;

        return glm::normalize(glm::vec3(ax * n.x, n.y, ay * n.z));
    }

    float GgxVndfAnisotropicPdf(const glm::vec3 &wi, const glm::vec3 &wm, const glm::vec3 &wo, float ax, float ay) {
        float absDotNL = std::abs(std::cos(wi.y));
        float absDotLH = std::abs(glm::dot(wm, wi));

        float G1 = Bsdf::SeparableSmithGGXG1(wo, wm, ax, ay);
        float D = Bsdf::GgxAnisotropicD(wm, ax, ay);

        return G1 * absDotLH * D / absDotNL;
    }

    void GgxVndfAnisotropicPdf(const glm::vec3 &wi, const glm::vec3 &wm, const glm::vec3 &wo, float ax, float ay,
                               float &forwardPdfW, float &reversePdfW) {
        float D = Bsdf::GgxAnisotropicD(wm, ax, ay);

        float absDotNL = std::abs(std::cos(wi.y));
        float absDotHL = std::abs(glm::dot(wm, wi));
        float G1v = Bsdf::SeparableSmithGGXG1(wo, wm, ax, ay);
        forwardPdfW = G1v * absDotHL * D / absDotNL;

        float absDotNV = std::abs(std::cos(wo.y));
        float absDotHV = std::abs(glm::dot(wm, wo));
        float G1l = Bsdf::SeparableSmithGGXG1(wi, wm, ax, ay);
        reversePdfW = G1l * absDotHV * D / absDotNV;
    }
}