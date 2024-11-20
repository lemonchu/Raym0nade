#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <glm/glm.hpp>
#include "material.h"

using vec2 = glm::vec<2, float>;
using vec3 = glm::vec<3, float>;
using vec4 = glm::vec<4, float>;

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

struct VertexData {
    vec2 uv;
    vec3 normal;

    VertexData();
    VertexData(const vec2 &uv, const vec3 &normal);
};

struct Face {
    vec3 v[3];
    VertexData* data[3];
    Material *texture;

    [[nodiscard]] vec3 center() const;
    [[nodiscard]] Box aabb() const;
};

struct HitRecord {
    float t_min, t_max;
    const Face* face;
    HitRecord();
    HitRecord(float t_min, float t_max);
};

bool RayTriangleIntersection(const Ray& ray, const Face& face, HitRecord& hit);

const float eps_edge = 1e-3f;
const float eps_zero = 1e-6f;

#endif // GEOMETRY_H