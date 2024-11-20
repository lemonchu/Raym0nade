#include <algorithm>
#include <iostream>
#include "kdt.h"

KDT_Node::KDT_Node() {
    son[0] = son[1] = nullptr;
}

const unsigned int LeafBagSize = 16;

template<int axis>
bool cmp(const Face &A, const Face &B) {
    return A.center()[axis] < B.center()[axis];
}

bool (*cmpFunc[3])(const Face &, const Face &) = {cmp<0>, cmp<1>, cmp<2>};

KDT::KDT() : buffer(nullptr), root(nullptr), cur(0) {}

KDT_Node *KDT::newNode() { return buffer + (cur++); }

void KDT::dfs_build(KDT_Node *&u, Face *faceL, Face *faceR) {
    u = newNode();
    if (faceR - faceL <= LeafBagSize) {
        u->faceL = faceL;
        u->faceR = faceR;
        u->box = Box(vec3(INFINITY),vec3(-INFINITY));
        for (Face *ptr = u->faceL; ptr != u->faceR; ++ptr)
            u->box = u->box + ptr->aabb();
        return;
    }
    vec3 Em = vec3(0), Em2 = vec3(0);
    for (Face *ptr = faceL; ptr != faceR; ++ptr) {
        vec3 m = (ptr->v[0] + ptr->v[1] + ptr->v[2]);
        Em += m;
        Em2 += m * m;
    }
    vec3 D = Em2 - Em * Em / (float)(faceR - faceL);
    int axis = 0;
    if (D[1] > D[0]) axis = 1;
    if (D[2] > D[0] && D[2] > D[1]) axis = 2;
    Face *faceM = faceL + (faceR - faceL) / 2;
    std::nth_element(faceL, faceM, faceR, cmpFunc[axis]);
    dfs_build(u->son[0], faceL, faceM);
    dfs_build(u->son[1], faceM, faceR);
    u->box = u->son[0]->box + u->son[1]->box;
}

void KDT::build(std::vector<Face> &triangles) {
    buffer = new KDT_Node[triangles.size() / LeafBagSize * 4];
    dfs_build(root, &triangles.front(), &triangles.back() + 1);
    std::cout << "build(cur): " << cur << std::endl;
}


void KDT::dfs_rayHit(KDT_Node *u) {
    if (u->son[0] == nullptr) {
        for (Face *ptr = u->faceL; ptr != u->faceR; ++ptr)
            RayTriangleIntersection(ray, *ptr, closest_hit);
        return;
    }

    float
        tL0 = closest_hit.t_min,
        tR0 = closest_hit.t_max,
        tL1 = closest_hit.t_min,
        tR1 = closest_hit.t_max;
    rayInBox(ray, u->son[0]->box, tL0, tR0);
    rayInBox(ray, u->son[1]->box, tL1, tR1);

    if (tL0 < tL1) {
        if (tL0 < tR0)
            dfs_rayHit(u->son[0]);
        if (tL1 < tR1 && tL1 < closest_hit.t_max)
            dfs_rayHit(u->son[1]);
    } else {
        if (tL1 < tR1)
            dfs_rayHit(u->son[1]);
        if (tL0 < tR0 && tL0 < closest_hit.t_max)
            dfs_rayHit(u->son[0]);
    }
}

HitRecord KDT::rayHit(Ray _ray) {
    ray = _ray;
    closest_hit = HitRecord(eps_zero, INFINITY);
    dfs_rayHit(root);
    return closest_hit;
}

KDT::~KDT() {
    delete[] buffer;
}