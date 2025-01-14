#include <algorithm>
#include <iostream>
#include "bvh.h"

BVH_Node::BVH_Node() : faceL(0), faceR(0) {}

template<int axis>
bool cmp(const Face &A, const Face &B) {
    return A.center()[axis] < B.center()[axis];
}

bool (*cmpFunc[3])(const Face &, const Face &) = {cmp<0>, cmp<1>, cmp<2>};

BVH::BVH() : node(nullptr), faces(nullptr) {}

const unsigned int LeafBagSize = 10;

void BVH::dfs_build(int u, int faceL, int faceR) {
    if (faceR - faceL <= LeafBagSize) {
        node[u].faceL = faceL;
        node[u].faceR = faceR;
        node[u].box = Box(vec3(INFINITY),vec3(-INFINITY));
        for (int i = faceL; i < faceR; i++)
            node[u].box = node[u].box + faces[i].aabb();
        return;
    }
    vec3 Em = vec3(0), Em2 = vec3(0);
    for (int i = faceL; i < faceR; i++) {
        vec3 m = faces[i].center();
        Em += m;
        Em2 += m * m;
    }
    vec3 D = Em2 - Em * Em / (float)(faceR - faceL);
    int axis = 0;
    if (D[1] > D[0]) axis = 1;
    if (D[2] > D[0] && D[2] > D[1]) axis = 2;
    int faceM = (faceR + faceL) / 2;
    std::nth_element(faces+faceL, faces+faceM, faces+faceR, cmpFunc[axis]);
    dfs_build(u<<1, faceL, faceM);
    dfs_build(u<<1|1, faceM, faceR);
    node[u].box = node[u<<1].box + node[u<<1|1].box;
}

int nodeCount(int u, int n) {
    return (n <= LeafBagSize) ? u : nodeCount(u<<1|1, (n+1)>>1);
}

void BVH::build(std::vector<Face> &faces0) {
    int size = nodeCount(1, int(faces0.size())) + 1;
    node = new BVH_Node[size];
    faces = &faces0.front();
    dfs_build(1, 0, int(faces0.size()));
    std::cout << "BVH has built with size " << size << std::endl;
}

void BVH::dfs_rayHit(int u, const Ray &ray, HitRecord &closest_hit) const {
    if (node[u].faceR) {
        for (int i = node[u].faceL; i < node[u].faceR; i++) {
            const Face &face = faces[i];
            float t = RayTriangleIntersection(ray, face.v[0], face.v[1], face.v[2]);
            if (closest_hit.t_min < t && t < closest_hit.t_max) {
                closest_hit.t_max = t;
                closest_hit.face = &face;
            }
        }
        return;
    }

    float
        tL0 = closest_hit.t_min,
        tR0 = closest_hit.t_max,
        tL1 = closest_hit.t_min,
        tR1 = closest_hit.t_max;
    rayInBox(ray, node[u<<1].box, tL0, tR0);
    rayInBox(ray, node[u<<1|1].box, tL1, tR1);

    if (tL0 < tL1) {
        if (tL0 < tR0)
            dfs_rayHit(u<<1, ray, closest_hit);
        if (tL1 < tR1 && tL1 < closest_hit.t_max)
            dfs_rayHit(u<<1|1, ray, closest_hit);
    } else {
        if (tL1 < tR1)
            dfs_rayHit(u<<1|1, ray, closest_hit);
        if (tL0 < tR0 && tL0 < closest_hit.t_max)
            dfs_rayHit(u<<1, ray, closest_hit);
    }
}

void BVH::rayHit(const Ray &ray, HitRecord &closest_hit) const {
    dfs_rayHit(1, ray, closest_hit);
}

BVH::~BVH() {
    delete[] node;
}