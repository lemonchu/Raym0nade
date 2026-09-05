#pragma once

#include <cstddef>
#include <cstdint>
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

class DirectLightSamplingScratch {
    // One instance may serve one concurrent sampling stream. Model identities keep cache
    // invalidation safe even when object storage addresses are reused.
private:
    friend void sampleDirectLight(
        const Bsdf& bsdf,
        const Model& model,
        Generator& generator,
        int sampleCount,
        std::vector<LightSample>& samples,
        DirectLightSamplingScratch& scratch);

    std::vector<double> areaLightWeights_;
    std::vector<double> areaLightCumulativeWeights_;
    std::uint64_t areaLightModelIdentity_{0U};
    vec3 areaLightPosition_{0.0F};
    std::size_t areaLightCount_{0U};
    double areaLightTotalWeight_{0.0};
    int areaLightLastValidIndex_{-1};
    bool areaLightWeightsInitialized_{false};
};

inline constexpr float kMinimumLightDistanceSquared = 1.0e-6F;

namespace sampling_detail {

// Exposed so exact boundary and zero-width-bin semantics remain regression tested.
[[nodiscard]] std::size_t firstCumulativeWeightAbove(
    const std::vector<double>& cumulativeWeights, double target) noexcept;

}  // namespace sampling_detail

void sampleDirectLight(
    const Bsdf& bsdf,
    const Model& model,
    Generator& generator,
    int sampleCount,
    std::vector<LightSample>& samples);

void sampleDirectLight(
    const Bsdf& bsdf,
    const Model& model,
    Generator& generator,
    int sampleCount,
    std::vector<LightSample>& samples,
    DirectLightSamplingScratch& scratch);

}  // namespace raym0nade
