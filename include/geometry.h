#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <glm/glm.hpp>
#include "texture.h"

using vec3 = glm::vec<3, float>;
using vec2 = glm::vec<2, float>;

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

struct Triangle {
    vec3 v[3];
    vec2 uv[3];
    Texture *texture;

    Box aabb();
};

struct HitRecord {
    float t_min, t_max;
    Triangle* tri;
    HitRecord();
    HitRecord(float t_min, float t_max);
};

bool RayTriangleIntersection(const Ray& ray, Triangle& tri, HitRecord& hit);

const float eps_edge = 1e-3f;
const float eps_zero = 1e-6f;

#endif // GEOMETRY_H