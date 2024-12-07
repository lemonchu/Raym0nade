#include "geometry.h"

bool isnan(vec3 v) {
    return std::isnan(v.x) || std::isnan(v.y) || std::isnan(v.z);
}


Box::Box() = default;

Box::Box(const vec3 &v0, const vec3 &v1) : v0(v0), v1(v1) {}

Box operator + (const Box &A, const Box &B) {
    return Box(
            vec3(std::fmin(A.v0[0], B.v0[0]), std::fmin(A.v0[1], B.v0[1]), std::fmin(A.v0[2], B.v0[2])),
            vec3(std::fmax(A.v1[0], B.v1[0]), std::fmax(A.v1[1], B.v1[1]), std::fmax(A.v1[2], B.v1[2]))
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
                tL = std::fmax(tL, (box.v0[i] - ray.origin[i]) * invD);
                tR = std::fmin(tR, (box.v1[i] - ray.origin[i]) * invD);
            } else {
                tL = std::fmax(tL, (box.v1[i] - ray.origin[i]) * invD);
                tR = std::fmin(tR, (box.v0[i] - ray.origin[i]) * invD);
            }
            tR += eps_zero;
            if (tL > tR)
                return ;
        }
    }
}

float RayTriangleIntersection(const Ray& ray, const vec3 &v0, const vec3 &v1, const vec3 &v2) {
    // Möller-Trumbore intersection algorithm
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;
    vec3 h = cross(ray.direction, edge2);
    float a = dot(edge1, h);

    if (std::abs(a) / length(edge1) < eps_zero)
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

vec3 barycentric(const vec3& v0, const vec3& v1, const vec3& v2, const vec3& P) {
    vec3 v0P = P - v0, v01 = v1 - v0, v02 = v2 - v0;
    float area = length(cross(v01, v02));
    float alpha = length(cross(v0P, v02)) / area;
    float beta = length(cross(v01, v0P)) / area;
    return vec3(1.0f - alpha - beta, alpha, beta);
}

void getTangentSpace(const vec3 &normal, vec3 &tangent, vec3 &bitangent) {
    vec3 v0 =
            abs(normal.x) < 0.8f ?
            vec3(1.0f, 0.0f, 0.0f) :
            vec3(0.0f, 1.0f, 0.0f);
    tangent = normalize(cross(v0, normal));
    bitangent = cross(normal, tangent);
}

const float eps_direction = 1e-3f;

void getTangentSpaceWithInDir(const vec3 &normal, const vec3 &inDir, vec3 &tangent, vec3 &bitangent) {
    float cosTheta = dot(normal, inDir);
    if (cosTheta < - 1.0f + eps_direction) {
        getTangentSpace(normal, tangent, bitangent);
        return ;
    }
    bitangent = normalize(cross(normal, inDir));
    tangent = cross(bitangent, normal);
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