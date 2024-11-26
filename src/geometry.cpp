#include "geometry.h"

Box::Box() {}

Box::Box(const vec3 &v0, const vec3 &v1) : v0(v0), v1(v1) {}

Box operator + (const Box &A, const Box &B) {
    return Box(
            vec3(std::min(A.v0[0], B.v0[0]), std::min(A.v0[1], B.v0[1]), std::min(A.v0[2], B.v0[2])),
            vec3(std::max(A.v1[0], B.v1[0]), std::max(A.v1[1], B.v1[1]), std::max(A.v1[2], B.v1[2]))
    );
}

void rayInBox(const Ray &ray, const Box &box, float &tL, float &tR) {
    for (int i = 0; i < 3; ++i) {
        if (std::abs(ray.direction[i]) < eps_zero) {
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
            tR += eps_zero;
            if (tL > tR)
                return ;
        }
    }
}

float RayTriangleIntersection(const Ray& ray, const vec3 &v0, const vec3 &v1, const vec3 &v2) {
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;
    vec3 h = cross(ray.direction, edge2);
    float a = dot(edge1, h);

    if (a <= 0.0f)
        return INFINITY;

    float f = 1.0f / a;
    vec3 s = ray.origin - v0;
    float u = f * dot(s, h);

    if (u < 0.0f || u > 1.0f)
        return INFINITY;

    vec3 q = cross(s, edge1);
    float v = f * dot(ray.direction, q);

    if (v < 0.0f || u + v > 1.0f)
        return INFINITY;

    return f * dot(edge2, q);
}

vec3 barycentric(const vec3& A, const vec3& B, const vec3& C, const vec3& P) {
    vec3 v0 = B - A;
    vec3 v1 = C - A;
    vec3 v2 = P - A;
    float d00 = dot(v0, v0);
    float d01 = dot(v0, v1);
    float d11 = dot(v1, v1);
    float d20 = dot(v2, v0);
    float d21 = dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;
    return vec3(u, v, w);
}

#include <iostream>

void getTangentSpace(const vec3 &normal, vec3 &tangent, vec3 &bitangent) {
    vec3 v0 =
            abs(normal.x) < 0.8f ?
            vec3(1.0f, 0.0f, 0.0f) :
            vec3(0.0f, 1.0f, 0.0f);
    tangent = normalize(cross(v0, normal));
    bitangent = cross(normal, tangent);
}

void tangentTransform(const vec3 &normal, vec3 &v) {
    vec3 v0 =
            abs(normal.x) < 0.8f ?
            vec3(1.0f, 0.0f, 0.0f) :
            vec3(0.0f, 1.0f, 0.0f);
    vec3 tangent, bitangent;
    getTangentSpace(normal, tangent, bitangent);
    v = v[0] * tangent + v[1] * bitangent + v[2] * normal;
}