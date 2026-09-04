#pragma once

#include <cstddef>
#include <vector>

#include "raym0nade/component.hpp"

namespace raym0nade {

class Bvh {
public:
    Bvh() = default;
    Bvh(const Bvh&) = delete;
    Bvh& operator=(const Bvh&) = delete;
    Bvh(Bvh&&) = delete;
    Bvh& operator=(Bvh&&) = delete;

    void build(std::vector<Face>& faces);
    void intersect(const Ray& ray, HitRecord& closestHit) const noexcept;

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::size_t nodeCount() const noexcept;

private:
    struct Node {
        Box bounds;
        std::size_t firstFace{0};
        std::size_t faceCount{0};
    };

    void buildNode(std::size_t nodeIndex, std::size_t firstFace, std::size_t faceEnd);
    void intersectNode(
        std::size_t nodeIndex,
        const Ray& ray,
        float nodeMinimum,
        float nodeMaximum,
        HitRecord& closestHit) const noexcept;

    std::vector<Node> nodes_;
    Face* faces_{nullptr};
};

}  // namespace raym0nade
