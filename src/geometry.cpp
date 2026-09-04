#include "raym0nade/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace raym0nade {
namespace {

constexpr float kDirectionEpsilon = 1.0e-7F;

vec3 normalizedOrZero(const vec3& value) noexcept {
    if (!isFinite(value)) {
        return vec3{0.0F};
    }

    const float scale = std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
    if (!(scale > std::numeric_limits<float>::min())) {
        return vec3{0.0F};
    }

    const vec3 scaled = value / scale;
    const float lengthSquared = glm::dot(scaled, scaled);
    if (!(lengthSquared > 0.0F) || !isFinite(lengthSquared)) {
        return vec3{0.0F};
    }

    return scaled / std::sqrt(lengthSquared);
}

}  // namespace

bool isFinite(float value) noexcept {
    return std::isfinite(value);
}

bool isFinite(const vec2& value) noexcept {
    return isFinite(value.x) && isFinite(value.y);
}

bool isFinite(const vec3& value) noexcept {
    return isFinite(value.x) && isFinite(value.y) && isFinite(value.z);
}

vec3 safeNormalize(const vec3& value, const vec3& fallback) noexcept {
    const vec3 normalized = normalizedOrZero(value);
    if (glm::dot(normalized, normalized) > 0.0F) {
        return normalized;
    }
    return normalizedOrZero(fallback);
}

float safeSqrt(float value) noexcept {
    return value > 0.0F ? std::sqrt(value) : 0.0F;
}

vec3 positivePow(vec3 value, float exponent) noexcept {
    if (!isFinite(exponent) || exponent < 0.0F) {
        return vec3{0.0F};
    }
    value = glm::max(value, vec3{0.0F});
    const vec3 result{
        std::pow(value.x, exponent),
        std::pow(value.y, exponent),
        std::pow(value.z, exponent),
    };
    return isFinite(result) ? result : vec3{0.0F};
}

Box::Box(const vec3& minimumValue, const vec3& maximumValue)
    : minimum(minimumValue), maximum(maximumValue) {}

Box merge(const Box& lhs, const Box& rhs) noexcept {
    return Box{
        vec3{
            std::fmin(lhs.minimum.x, rhs.minimum.x),
            std::fmin(lhs.minimum.y, rhs.minimum.y),
            std::fmin(lhs.minimum.z, rhs.minimum.z),
        },
        vec3{
            std::fmax(lhs.maximum.x, rhs.maximum.x),
            std::fmax(lhs.maximum.y, rhs.maximum.y),
            std::fmax(lhs.maximum.z, rhs.maximum.z),
        },
    };
}

bool intersect(const Ray& ray, const Box& box, float& tMinimum, float& tMaximum) noexcept {
    if (!isFinite(ray.origin) || !isFinite(ray.direction) || !isFinite(box.minimum) ||
        !isFinite(box.maximum) || std::isnan(tMinimum) || std::isnan(tMaximum) ||
        tMinimum > tMaximum) {
        return false;
    }

    double nearDistance = static_cast<double>(tMinimum);
    double farDistance = static_cast<double>(tMaximum);

    for (int axis = 0; axis < 3; ++axis) {
        const double minimumValue = static_cast<double>(box.minimum[axis]);
        const double maximumValue = static_cast<double>(box.maximum[axis]);
        if (minimumValue > maximumValue) {
            return false;
        }

        const double origin = static_cast<double>(ray.origin[axis]);
        const double direction = static_cast<double>(ray.direction[axis]);
        if (direction == 0.0) {
            if (origin < minimumValue || origin > maximumValue) {
                return false;
            }
            continue;
        }

        double axisNear = (minimumValue - origin) / direction;
        double axisFar = (maximumValue - origin) / direction;
        if (axisNear > axisFar) {
            std::swap(axisNear, axisFar);
        }

        nearDistance = std::max(nearDistance, axisNear);
        farDistance = std::min(farDistance, axisFar);
        if (nearDistance > farDistance) {
            return false;
        }
    }

    tMinimum = static_cast<float>(nearDistance);
    tMaximum = static_cast<float>(farDistance);
    return true;
}

float intersectTriangle(const Ray& ray, const vec3& v0, const vec3& v1, const vec3& v2) noexcept {
    constexpr float kDeterminantTolerance = 8.0F * std::numeric_limits<float>::epsilon();
    const float infinity = std::numeric_limits<float>::infinity();

    if (!isFinite(ray.origin) || !isFinite(ray.direction) || !isFinite(v0) || !isFinite(v1) ||
        !isFinite(v2)) {
        return infinity;
    }

    const vec3 edge1 = v1 - v0;
    const vec3 edge2 = v2 - v0;
    const vec3 normal = glm::cross(edge1, edge2);
    const float normalLengthSquared = glm::dot(normal, normal);
    if (!(normalLengthSquared > std::numeric_limits<float>::min()) ||
        !isFinite(normalLengthSquared)) {
        return infinity;
    }

    const vec3 perpendicular = glm::cross(ray.direction, edge2);
    const float determinant = glm::dot(edge1, perpendicular);
    const float determinantScale = std::sqrt(glm::dot(edge1, edge1) * glm::dot(perpendicular, perpendicular));
    if (!isFinite(determinantScale) || !(determinantScale > 0.0F) ||
        std::abs(determinant) <= kDeterminantTolerance * determinantScale) {
        return infinity;
    }

    const float inverseDeterminant = 1.0F / determinant;
    const vec3 originOffset = ray.origin - v0;
    const float u = inverseDeterminant * glm::dot(originOffset, perpendicular);
    if (u < 0.0F || u > 1.0F) {
        return infinity;
    }

    const vec3 crossOffset = glm::cross(originOffset, edge1);
    const float v = inverseDeterminant * glm::dot(ray.direction, crossOffset);
    if (v < 0.0F || u + v > 1.0F) {
        return infinity;
    }

    const float distance = inverseDeterminant * glm::dot(edge2, crossOffset);
    return isFinite(distance) ? distance : infinity;
}

vec3 barycentric(const vec3& v0, const vec3& v1, const vec3& v2, const vec3& point) noexcept {
    if (!isFinite(v0) || !isFinite(v1) || !isFinite(v2) || !isFinite(point)) {
        return vec3{0.0F};
    }

    const vec3 edge1 = v1 - v0;
    const vec3 edge2 = v2 - v0;
    const vec3 normal = glm::cross(edge1, edge2);
    const float denominator = glm::dot(normal, normal);
    if (!(denominator > std::numeric_limits<float>::min()) || !isFinite(denominator)) {
        return vec3{0.0F};
    }

    const vec3 offset = point - v0;
    const float second = glm::dot(glm::cross(offset, edge2), normal) / denominator;
    const float third = glm::dot(glm::cross(edge1, offset), normal) / denominator;
    const vec3 result{1.0F - second - third, second, third};
    return isFinite(result) ? result : vec3{0.0F};
}

void makeTangentSpace(const vec3& normal, vec3& tangent, vec3& bitangent) noexcept {
    const vec3 unitNormal = safeNormalize(normal, vec3{0.0F, 0.0F, 1.0F});
    const vec3 helper = std::abs(unitNormal.x) < 0.8F ? vec3{1.0F, 0.0F, 0.0F}
                                                        : vec3{0.0F, 1.0F, 0.0F};
    tangent = safeNormalize(glm::cross(helper, unitNormal), vec3{0.0F, 1.0F, 0.0F});
    bitangent = safeNormalize(glm::cross(unitNormal, tangent), vec3{1.0F, 0.0F, 0.0F});
}

void makeTangentSpace(
    const vec3& normal, const vec3& incoming, vec3& tangent, vec3& bitangent) noexcept {
    const vec3 unitNormal = safeNormalize(normal, vec3{0.0F, 0.0F, 1.0F});
    const vec3 candidate = glm::cross(unitNormal, incoming);
    if (!isFinite(candidate) || glm::dot(candidate, candidate) <= kDirectionEpsilon * kDirectionEpsilon) {
        makeTangentSpace(unitNormal, tangent, bitangent);
        return;
    }

    bitangent = safeNormalize(candidate);
    tangent = safeNormalize(glm::cross(bitangent, unitNormal));
}

}  // namespace raym0nade
