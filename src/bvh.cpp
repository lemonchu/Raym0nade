#include "raym0nade/bvh.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace raym0nade {
namespace {

constexpr std::size_t kLeafFaceCount = 10;

float centroidCoordinate(const Face& face, int axis) noexcept {
    const float coordinate = face.center()[axis];
    return std::isfinite(coordinate) ? coordinate : 0.0F;
}

}  // namespace

void Bvh::build(std::vector<Face>& faces) {
    nodes_.clear();
    faces_ = nullptr;
    if (faces.empty()) {
        return;
    }

    faces_ = faces.data();
    nodes_.emplace_back();
    buildNode(0, 0, faces.size());
}

void Bvh::buildNode(std::size_t nodeIndex, std::size_t firstFace, std::size_t faceEnd) {
    Box bounds;
    for (std::size_t index = firstFace; index < faceEnd; ++index) {
        bounds = merge(bounds, faces_[index].bounds());
    }

    const std::size_t faceCount = faceEnd - firstFace;
    if (faceCount <= kLeafFaceCount) {
        nodes_[nodeIndex].bounds = bounds;
        nodes_[nodeIndex].firstFace = firstFace;
        nodes_[nodeIndex].faceCount = faceCount;
        return;
    }

    vec3 centroidMinimum{std::numeric_limits<float>::infinity()};
    vec3 centroidMaximum{-std::numeric_limits<float>::infinity()};
    for (std::size_t index = firstFace; index < faceEnd; ++index) {
        const vec3 centroid{
            centroidCoordinate(faces_[index], 0),
            centroidCoordinate(faces_[index], 1),
            centroidCoordinate(faces_[index], 2),
        };
        centroidMinimum = glm::min(centroidMinimum, centroid);
        centroidMaximum = glm::max(centroidMaximum, centroid);
    }

    const vec3 extent = centroidMaximum - centroidMinimum;
    int axis = 0;
    if (extent.y > extent.x) {
        axis = 1;
    }
    if (extent.z > extent[axis]) {
        axis = 2;
    }

    const std::size_t middleFace = firstFace + faceCount / 2;
    std::nth_element(
        faces_ + firstFace,
        faces_ + middleFace,
        faces_ + faceEnd,
        [axis](const Face& lhs, const Face& rhs) {
            return centroidCoordinate(lhs, axis) < centroidCoordinate(rhs, axis);
        });

    const std::size_t firstChild = nodes_.size();
    nodes_.resize(nodes_.size() + 2);
    nodes_[nodeIndex].bounds = bounds;
    nodes_[nodeIndex].firstFace = firstChild;
    nodes_[nodeIndex].faceCount = 0;

    buildNode(firstChild, firstFace, middleFace);
    buildNode(firstChild + 1, middleFace, faceEnd);
}

void Bvh::intersectNode(
    std::size_t nodeIndex,
    const Ray& ray,
    float nodeMinimum,
    float nodeMaximum,
    HitRecord& closestHit) const noexcept {
    if (nodeIndex >= nodes_.size() || nodeMinimum >= closestHit.tMaximum) {
        return;
    }

    const Node& node = nodes_[nodeIndex];
    if (node.faceCount > 0) {
        const std::size_t faceEnd = node.firstFace + node.faceCount;
        for (std::size_t index = node.firstFace; index < faceEnd; ++index) {
            const Face& face = faces_[index];
            const float distance = intersectTriangle(
                ray, face.vertices[0], face.vertices[1], face.vertices[2]);
            if (distance > closestHit.tMinimum && distance < closestHit.tMaximum) {
                closestHit.tMaximum = distance;
                closestHit.face = &face;
            }
        }
        return;
    }

    const std::size_t leftChild = node.firstFace;
    const std::size_t rightChild = leftChild + 1;
    if (rightChild >= nodes_.size()) {
        return;
    }

    const float traversalMinimum = std::max(nodeMinimum, closestHit.tMinimum);
    const float traversalMaximum = std::min(nodeMaximum, closestHit.tMaximum);
    float leftMinimum = traversalMinimum;
    float leftMaximum = traversalMaximum;
    float rightMinimum = traversalMinimum;
    float rightMaximum = traversalMaximum;
    const bool intersectsLeft = raym0nade::intersect(
        ray, nodes_[leftChild].bounds, leftMinimum, leftMaximum);
    const bool intersectsRight = raym0nade::intersect(
        ray, nodes_[rightChild].bounds, rightMinimum, rightMaximum);

    if (intersectsLeft && intersectsRight) {
        if (leftMinimum <= rightMinimum) {
            intersectNode(leftChild, ray, leftMinimum, leftMaximum, closestHit);
            if (rightMinimum < closestHit.tMaximum) {
                intersectNode(
                    rightChild,
                    ray,
                    rightMinimum,
                    std::min(rightMaximum, closestHit.tMaximum),
                    closestHit);
            }
        } else {
            intersectNode(rightChild, ray, rightMinimum, rightMaximum, closestHit);
            if (leftMinimum < closestHit.tMaximum) {
                intersectNode(
                    leftChild,
                    ray,
                    leftMinimum,
                    std::min(leftMaximum, closestHit.tMaximum),
                    closestHit);
            }
        }
    } else if (intersectsLeft) {
        intersectNode(leftChild, ray, leftMinimum, leftMaximum, closestHit);
    } else if (intersectsRight) {
        intersectNode(rightChild, ray, rightMinimum, rightMaximum, closestHit);
    }
}

void Bvh::intersect(const Ray& ray, HitRecord& closestHit) const noexcept {
    if (empty()) {
        return;
    }

    float rootMinimum = closestHit.tMinimum;
    float rootMaximum = closestHit.tMaximum;
    if (raym0nade::intersect(ray, nodes_.front().bounds, rootMinimum, rootMaximum)) {
        intersectNode(0, ray, rootMinimum, rootMaximum, closestHit);
    }
}

bool Bvh::empty() const noexcept {
    return nodes_.empty() || faces_ == nullptr;
}

std::size_t Bvh::nodeCount() const noexcept {
    return nodes_.size();
}

}  // namespace raym0nade
