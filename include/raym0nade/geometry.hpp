#pragma once

#include <cmath>
#include <limits>

#include <glm/glm.hpp>

namespace raym0nade {

using vec2 = glm::vec2;
using vec3 = glm::vec3;
using vec4 = glm::vec4;

inline constexpr float kRayEpsilon = 1.0e-4F;
inline constexpr float kPi = 3.14159265358979323846F;
inline constexpr vec3 kLuminanceWeights{0.3F, 0.6F, 0.1F};

[[nodiscard]] bool isFinite(float value) noexcept;
[[nodiscard]] bool isFinite(const vec2& value) noexcept;
[[nodiscard]] bool isFinite(const vec3& value) noexcept;
[[nodiscard]] vec3 safeNormalize(const vec3& value, const vec3& fallback = vec3{0.0F}) noexcept;
[[nodiscard]] float safeSqrt(float value) noexcept;
[[nodiscard]] vec3 positivePow(vec3 value, float exponent) noexcept;

struct Ray {
    vec3 origin{0.0F};
    vec3 direction{0.0F, 0.0F, 1.0F};
};

struct RayDifferential {
    vec3 dPdx{0.0F};
    vec3 dPdy{0.0F};
    vec3 dDdx{0.0F};
    vec3 dDdy{0.0F};
};

struct Box {
    vec3 minimum{std::numeric_limits<float>::infinity()};
    vec3 maximum{-std::numeric_limits<float>::infinity()};

    Box() = default;
    Box(const vec3& minimum, const vec3& maximum);
};

[[nodiscard]] Box merge(const Box& lhs, const Box& rhs) noexcept;
[[nodiscard]] bool intersect(const Ray& ray, const Box& box, float& tMinimum, float& tMaximum) noexcept;
[[nodiscard]] float intersectTriangle(
    const Ray& ray, const vec3& v0, const vec3& v1, const vec3& v2) noexcept;
[[nodiscard]] vec3 barycentric(
    const vec3& v0, const vec3& v1, const vec3& v2, const vec3& point) noexcept;

void makeTangentSpace(const vec3& normal, vec3& tangent, vec3& bitangent) noexcept;
void makeTangentSpace(
    const vec3& normal, const vec3& incoming, vec3& tangent, vec3& bitangent) noexcept;

}  // namespace raym0nade
