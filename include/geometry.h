#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <glm/glm.hpp>
#include <cmath>

#define RAY_DEBUG

using vec2 = glm::vec<2, float>;
using vec3 = glm::vec<3, float>;
using vec4 = glm::vec<4, float>;

const float eps_zero = 1e-4f;

const vec3 RGB_Weight = vec3(0.3f, 0.6f, 0.1f);

bool finite(vec2 v);
bool finite(vec3 v);

float sqrt_s(float x);
vec3 pow_s(vec3 v, float p);

struct Ray {
    vec3 origin, direction;
};

struct RayDifferential {
    vec3 dPdx = glm::vec3{0.0f}, dPdy = glm::vec3{0.0f};
    vec3 dDdx = glm::vec3{0.0f}, dDdy = glm::vec3{0.0f};
};

struct Box {
    vec3 v0, v1;
    Box();
    Box(const vec3 &v0, const vec3 &v1);
};

Box operator + (const Box &A, const Box &B);

void rayInBox(const Ray &ray, const Box &box, float &tL, float &tR);

float RayTriangleIntersection(const Ray& ray, const vec3 &v0, const vec3 &v1, const vec3 &v2);

glm::vec3 barycentric(const glm::vec3& A, const glm::vec3& B, const glm::vec3& C, const glm::vec3& P);

void getTangentSpace(const vec3 &normal, vec3 &tangent, vec3 &bitangent);
void getTangentSpaceWithInDir(const vec3 &normal, const vec3 &inDir, vec3 &tangent, vec3 &bitangent);

#endif // GEOMETRY_H