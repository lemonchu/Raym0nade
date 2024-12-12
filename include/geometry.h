#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <glm/glm.hpp>
#include <cmath>

using vec2 = glm::vec<2, float>;
using vec3 = glm::vec<3, float>;
using vec4 = glm::vec<4, float>;

#ifndef M_PI
const float M_PI = 3.14159265359f;
#endif

bool isnan(vec3 v);

struct Ray {
    vec3 origin, direction;
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
void tangentTransform(const vec3 &normal, vec3 &v);

const float eps_zero = 1e-4f;

#endif // GEOMETRY_H