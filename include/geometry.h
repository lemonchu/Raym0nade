#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <glm/glm.hpp>
#include <iostream>
#include "texture.h"

using vec3 = glm::vec<3, float>;
using vec2 = glm::vec<2, float>;

struct Ray {
    vec3 origin, direction;
};

struct Box {
    vec3 v0, v1;
    Box() {}
    Box(const vec3 &v0, const vec3 &v1) : v0(v0), v1(v1) {}
};

Box operator + (const Box &A, const Box &B) {
    return Box(
        vec3(std::min(A.v0[0], B.v0[0]), std::min(A.v0[1], B.v0[1]), std::min(A.v0[2], B.v0[2])),
        vec3(std::max(A.v1[0], B.v1[0]), std::max(A.v1[1], B.v1[1]), std::max(A.v1[2], B.v1[2]))
    );
}

void rayInBox(const Ray &ray, const Box &box, float &tL, float &tR) {

    const float eps = 1e-6f;

    for (int i = 0; i < 3; ++i) {
        if (std::abs(ray.direction[i]) < eps) {
            if (ray.origin[i] < box.v0[i] || ray.origin[i] > box.v1[i]) {
                tR = -1.0;
                return ;
            }
        } else {
            float invD = 1.0f / ray.direction[i];
            if (invD >= 0) {
                tL = std::max(tL, (box.v0[i] - ray.origin[i]) * invD);
                tR = std::min(tR, (box.v1[i] - ray.origin[i]) * invD);
            } else {
                tL = std::max(tL, (box.v1[i] - ray.origin[i]) * invD);
                tR = std::min(tR, (box.v0[i] - ray.origin[i]) * invD);
            }
            if (tL > tR)
                return ;
        }
    }
}

struct Triangle {
    vec3 v[3];
    vec2 uv[3];
    Texture* texture;
    Box aabb() {
        return Box(
            vec3(
                std::min(v[0][0],std::min(v[1][0],v[2][0])),
                std::min(v[0][1],std::min(v[1][1],v[2][1])),
                std::min(v[0][2],std::min(v[1][2],v[2][2]))
            ),
            vec3(
                std::max(v[0][0],std::max(v[1][0],v[2][0])),
                std::max(v[0][1],std::max(v[1][1],v[2][1])),
                std::max(v[0][2],std::max(v[1][2],v[2][2]))
            )
        );
    }
};

const float eps = 1e-4;

struct HitRecord {
    float t_min, t_max;
    Triangle* tri;
    HitRecord(float t_min, float t_max) : t_min(t_min), t_max(t_max), tri(nullptr) {}
    HitRecord() : t_min(eps), t_max(INFINITY), tri(nullptr) {}
};

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
        hit.t_max = t;
        hit.tri = &tri;
        return true;
    }

    return false;
}

#endif // GEOMETRY_H