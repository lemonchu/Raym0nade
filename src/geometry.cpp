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
            tR += eps_edge;
            if (tL > tR)
                return ;
        }
    }
}

VertexData::VertexData() {}

VertexData::VertexData(const vec2 &uv, const vec3 &normal, const vec4 &color) : uv(uv), normal(normal), color(color) {}

vec3 Face::center() const {
    return (v[0] + v[1] + v[2]) / 3.0f;
}

Box Face::aabb() const {
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

HitRecord::HitRecord() : t_min(eps_zero), t_max(INFINITY), face(nullptr) {}

HitRecord::HitRecord(float t_min, float t_max) : t_min(t_min), t_max(t_max), face(nullptr) {}

bool RayTriangleIntersection(const Ray& ray, const Face& face, HitRecord& hit) {
    vec3
        edge1 = face.v[1] - face.v[0],
        edge2 = face.v[2] - face.v[0],
        h = cross(ray.direction, edge2);
    float a = dot(edge1, h);

    if (abs(a) < eps_zero)
        return false;

    float f = 1.0 / a;
    vec3 s = ray.origin - face.v[0];
    float u = f * dot(s, h);

    if (u < -eps_edge || u > 1.0 + eps_edge)
        return false;

    vec3 q = cross(s, edge1);
    float v = f * dot(ray.direction, q);

    if (v < -eps_edge || u + v > 1.0 + eps_edge)
        return false;

    float t = f * dot(edge2, q);

    if (t < hit.t_max && t > hit.t_min) {
        hit.t_max = t;
        hit.face = &face;
        return true;
    }

    return false;
}