#pragma once

#include <vector>

#include "raym0nade/model.hpp"

namespace raym0nade {

class Bsdf {
public:
    Bsdf(const vec3& incomingDirection, const vec3& hitPosition);
    Bsdf(const vec3& incomingDirection, const HitInfo& hitInfo);

    [[nodiscard]] vec3 evaluate(const vec3& outgoingDirection) const;
    void refractPrecisely(vec3& outgoingDirection, float& fresnel) const noexcept;
    void sampleReflection(Generator& generator, vec3& direction, vec3& throughput, int& failures) const;
    void sampleTransmission(Generator& generator, vec3& direction, vec3& throughput, int& failures) const;

    vec3 incomingDirection;
    HitInfo surface;

private:
    void sampleMicrofacet(
        Generator& generator, vec3& outgoingDirection, vec3& throughput, float& pdf, int& failures) const;
    void sampleCosine(
        Generator& generator, vec3& outgoingDirection, vec3& throughput, float& pdf, int& failures) const;

    [[nodiscard]] vec3 evaluateReflection(const vec3& outgoingDirection) const;
    [[nodiscard]] vec3 evaluateTransmission(const vec3& outgoingDirection) const;
};

struct LightSample {
    vec3 throughput{0.0F};
    vec3 radiance{0.0F};
    float weight{0.0F};
};

inline constexpr float kMinimumLightDistanceSquared = 1.0e-6F;

void sampleDirectLight(
    const Bsdf& bsdf,
    const Model& model,
    Generator& generator,
    int sampleCount,
    std::vector<LightSample>& samples);

}  // namespace raym0nade
