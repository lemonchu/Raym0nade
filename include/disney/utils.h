#ifndef RAYM0NADE_UTILS_H
#define RAYM0NADE_UTILS_H

#include <cmath>
#include <algorithm>
#include <glm/glm.hpp>
#include <cstdlib>

inline float ReLU(float x) {
    return x > 0.0f ? x : 0.0f;
}

inline float lerp(float a, float b, float t) {
    return a + t * (b - a);
}

inline glm::vec3 lerp(const glm::vec3 &a, const glm::vec3 &b, float t) {
    return a + t * (b - a);
}

enum SurfaceEventFlags
{
    eScatterEvent        = 0x01,
    eTransmissionEvent   = 0x02,
    eDiracEvent          = 0x04
};

enum MediumPhaseFunction
{
    eVacuum,
    eIsotropic
};

struct MediumParameters
{
    MediumPhaseFunction phaseFunction = eVacuum;
    glm::vec3 extinction = glm::vec3{0.0f};
};

struct BsdfSample {
    uint32_t flags;

    MediumParameters medium = MediumParameters();
    glm::vec3 reflectance = glm::vec3{0.0f};
    glm::vec3 wi = glm::vec3{0.0f};
    float forwardPdfW = 0.0f;
    float reversePdfW = 0.0f;
};

namespace Bsdf {
    float SeparableSmithGGXG1(const glm::vec3 &w, const glm::vec3 &wm, float ax, float ay);

    float SeparableSmithGGXG1(const glm::vec3 &w, float a);

    float HeightCorrelatedSmithGGXG2(const glm::vec3 &wo, const glm::vec3 &wi, float a);

    float GgxIsotropicD(const glm::vec3 &wm, float a);

    float GgxAnisotropicD(const glm::vec3 &wm, float ax, float ay);

    glm::vec3 SampleGgxVndf(glm::vec3 wo, float roughness, float u1, float u2);

    float GgxVndfPdf(const glm::vec3 &wo, const glm::vec3 &wm, const glm::vec3 &wi, float a);

    glm::vec3 SampleGgxVndfAnisotropic(const glm::vec3 &wo, float ax, float ay, float u1, float u2);

    float GgxVndfAnisotropicPdf(const glm::vec3 &wi, const glm::vec3 &wm, const glm::vec3 &wo, float ax, float ay);

    void GgxVndfAnisotropicPdf(const glm::vec3 &wi, const glm::vec3 &wm, const glm::vec3 &wo, float ax, float ay,
                               float &forwardPdfW, float &reversePdfW);
}

inline float CosTheta(const glm::vec3 &v) {
    return v.y;
}

inline float AbsCosTheta(const glm::vec3 &v) {
    return std::abs(v.y);
}

inline float Dot(const glm::vec3 &a, const glm::vec3 &b) {
    return glm::dot(a, b);
}

inline float Absf(float x) {
    return std::abs(x);
}

inline float AbsDot(const glm::vec3 &a, const glm::vec3 &b) {
    return std::abs(glm::dot(a, b));
}

inline glm::vec3 Reflect(const glm::vec3 &v, const glm::vec3 &n) {
    return glm::reflect(v, n);
}

inline float Square(float x) {
    return x * x;
}

inline glm::vec3 Sqrt(const glm::vec3 &v) {
    return glm::sqrt(v);
}

inline glm::vec3 MatrixMultiply(const glm::vec3 &v, const glm::mat3 &m) {
    return m * v;
}

inline glm::mat3 MatrixTranspose(const glm::mat3 &m) {
    return glm::transpose(m);
}

inline float Sign(float x) {
    return x < 0.0f ? -1.0f : 1.0f;
}

inline bool Transmit(glm::vec3 wm, glm::vec3 wi, float n, glm::vec3& wo){
    float c = Dot(wi, wm);
    if(c < 0.0f) {
        c = -c;
        wm = -wm;
    }

    float root = 1.0f - n * n * (1.0f - c * c);
    if(root <= 0)
        return false;

    wo = (n * c - std::sqrtf(root)) * wm - n * wi;
    return true;
}

#endif //RAYM0NADE_UTILS_H