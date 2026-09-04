#pragma once

#include <cstdint>
#include <filesystem>
#include <limits>
#include <random>
#include <vector>

#include "raym0nade/geometry.hpp"

namespace raym0nade {

class Material;

struct VertexData {
    vec2 uv{0.0F};
    vec3 normal{0.0F};
};

struct Face {
    vec3 vertices[3]{};
    const VertexData* vertexData[3]{};
    const Material* material{nullptr};

    [[nodiscard]] vec3 center() const noexcept;
    [[nodiscard]] Box bounds() const noexcept;
};

struct HitRecord {
    float tMinimum{kRayEpsilon};
    float tMaximum{std::numeric_limits<float>::infinity()};
    const Face* face{nullptr};

    HitRecord() = default;
    HitRecord(float tMinimum, float tMaximum);
};

class Generator {
public:
    explicit Generator(std::uint32_t seed = 0);

    [[nodiscard]] float operator()() noexcept;

private:
    std::mt19937 engine_;
};

class RandomDistribution {
public:
    void initialize(const std::vector<float>& weights);

    [[nodiscard]] int sample(Generator& generator) const noexcept;
    [[nodiscard]] int operator()(Generator& generator) const noexcept;
    [[nodiscard]] float pdf(std::size_t index) const noexcept;
    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;

private:
    std::vector<double> cumulativeWeights_;
    double totalWeight_{0.0};
};

struct LightObject {
    vec3 center{0.0F};
    vec3 color{0.0F};
    float power{0.0F};
    std::vector<Face> faces;
    RandomDistribution faceDistribution;
};

class SkyBox {
public:
    void load(const std::filesystem::path& filename);

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] vec3 radiance(const vec3& direction) const noexcept;
    [[nodiscard]] bool sample(
        Generator& generator, int& pixelIndex, vec3& direction, vec3& weightedRadiance) const noexcept;

private:
    void initializeDistribution();

    int width_{0};
    int height_{0};
    std::vector<vec3> radiance_;
    std::vector<float> solidAngles_;
    RandomDistribution distribution_;
};

}  // namespace raym0nade
