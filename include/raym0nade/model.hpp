#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <vector>

#include "raym0nade/bvh.hpp"
#include "raym0nade/material.hpp"
#include "raym0nade/scene_data.hpp"

namespace raym0nade {

class ModelBuilder;

struct HitInfo {
    vec3 shapeNormal{std::numeric_limits<float>::quiet_NaN()};
    vec3 surfaceNormal{std::numeric_limits<float>::quiet_NaN()};
    vec3 emission{0.0F};
    vec3 baseColor{0.0F};
    vec3 position{std::numeric_limits<float>::quiet_NaN()};
    float specular{0.04F};
    float roughness{0.8F};
    float metallic{0.0F};
    float opacity{1.0F};
    float eta{1.0F};
    int materialId{0};
    bool entering{true};
};

void getHitNormals(
    const Face& face,
    const vec3& incomingDirection,
    const vec3& barycentricCoordinates,
    vec3& shapeNormal,
    vec3& surfaceNormal,
    bool& entering) noexcept;
void getHitMaterial(
    const Face& face,
    const vec3& barycentricCoordinates,
    const vec3& positionDx,
    const vec3& positionDy,
    HitInfo& hitInfo);

class Model {
public:
    Model(
        const std::filesystem::path& modelDirectory,
        const std::filesystem::path& modelFilename,
        const std::filesystem::path& skyFilename);
    Model(const Model&) = delete;
    Model& operator=(const Model&) = delete;
    Model(Model&&) = delete;
    Model& operator=(Model&&) = delete;

    [[nodiscard]] HitRecord intersect(const Ray& ray) const noexcept;
    [[nodiscard]] bool occluded(const Ray& ray, float maximumDistance) const noexcept;

    [[nodiscard]] const std::vector<LightObject>& lights() const noexcept;
    [[nodiscard]] const SkyBox& sky() const noexcept;
    [[nodiscard]] const std::filesystem::path& modelPath() const noexcept;
    // Process-unique for transient cache keys.
    [[nodiscard]] std::uint64_t instanceIdentity() const noexcept;
    [[nodiscard]] std::size_t faceCount() const noexcept;
    [[nodiscard]] PackedSceneData packScene() const;

private:
    friend class ModelBuilder;

    std::uint64_t instanceIdentity_{0U};
    std::vector<Material> materials_;
    std::vector<Face> faces_;
    std::vector<VertexData> vertexData_;
    std::vector<LightObject> lights_;
    Bvh bvh_;
    std::filesystem::path modelPath_;
    std::filesystem::path skyPath_;
    SkyBox sky_;
};

}  // namespace raym0nade
