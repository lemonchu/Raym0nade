#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <glm/glm.hpp>
#include <iostream>
#include "texture.h"

using vec3 = glm::vec<3, float>;
using vec2 = glm::vec<2, float>;

struct Triangle {
    vec3 v[3];
    vec2 uv[3];
    Texture* texture;
};

struct Ray {
    vec3 origin, direction;
};

struct HitRecord {
    vec3 normal;
    float t_min, t_max;
    Triangle* tri;
    HitRecord(vec3 normal, float t_min, float t_max) : normal(normal), t_min(t_min), t_max(t_max) {}
};

const float eps = 1e-4;

bool RayTriangleIntersection(const Ray& ray, Triangle& tri, HitRecord &hit) {
    vec3
        edge1 = tri.v[1] - tri.v[0],
        edge2 = tri.v[2] - tri.v[0],
        h = cross(ray.direction, edge2);
    float a = dot(edge1, h);

    if (abs(a) < eps)
        return false;

    float f = 1.0 / a;
    vec3 s = ray.origin - tri.v[0];
    float u = f * dot(s, h);

    if (u < 0.0 || u > 1.0)
        return false;

    vec3 q = cross(s, edge1);
    float v = f * dot(ray.direction, q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    float t = f * dot(edge2, q);

    if (t < hit.t_max && t > hit.t_min) {
        hit.normal = normalize(cross(edge1, edge2));
        hit.t_max = t;
        hit.tri = &tri;
        return true;
    }

    return false;
}

#endif // GEOMETRY_H