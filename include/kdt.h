#ifndef KDT_H
#define KDT_H

#include <algorithm>
#include "geometry.h"

struct KDT_Node {
    Box box;
    KDT_Node *son[2];
    Triangle *triL, *triR;

    KDT_Node() {
        son[0] = son[1] = nullptr;
    }
};

const unsigned int LeafBagSize = 16;

bool cmpX(const Triangle &A, const Triangle &B) {
    return (A.v[0] + A.v[1] + A.v[2])[0] < (B.v[0] + B.v[1] + B.v[2])[0];
}

bool cmpY(const Triangle &A, const Triangle &B) {
    return (A.v[0] + A.v[1] + A.v[2])[1] < (B.v[0] + B.v[1] + B.v[2])[1];
}

bool cmpZ(const Triangle &A, const Triangle &B) {
    return (A.v[0] + A.v[1] + A.v[2])[2] < (B.v[0] + B.v[1] + B.v[2])[2];
}

bool (*cmpFunc[3])(const Triangle &, const Triangle &) = {cmpX, cmpY, cmpZ};

class KDT {
    KDT_Node *buffer, *root;
    int cur;
public:
    KDT() {}

    KDT_Node *newNode() { return buffer + (cur++); }

    void dfs_build(KDT_Node *&u, Triangle *triL, Triangle *triR, int depth) {
        u = newNode();
        if (triR - triL <= LeafBagSize) {
            u->triL = triL;
            u->triR = triR;
            u->box = (Box){vec3(INFINITY),vec3(-INFINITY)};
            for (Triangle *ptr = u->triL; ptr != u->triR; ++ptr)
                u->box = u->box + ptr->aabb();
            return;
        }
        Triangle *triM = triL + (triR - triL) / 2;
        std::nth_element(triL, triM, triR, cmpFunc[depth%3]);
        dfs_build(u->son[0], triL, triM, depth + 1);
        dfs_build(u->son[1], triM, triR, depth + 1);
        u->box = u->son[0]->box + u->son[1]->box;
    }

    void build(std::vector<Triangle> &triangles) {
        buffer = new KDT_Node[triangles.size() / LeafBagSize * 4];
        cur = 0;
        dfs_build(root, &triangles.front(), &triangles.back() + 1, 0);
        std::cout << "build(cur): " << cur << std::endl;
    }

    HitRecord closest_hit;
    Ray ray;
    int cnt;

    void dfs_rayHit(KDT_Node *u) {
        if (u->son[0] == nullptr) {
            for (Triangle *ptr = u->triL; ptr != u->triR; ++ptr) {
                cnt++;
                RayTriangleIntersection(ray, *ptr, closest_hit);
            }
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
            if (tL0<tR0)
                dfs_rayHit(u->son[0]);
            if (tL1<tR1 && tL1<closest_hit.t_max)
                dfs_rayHit(u->son[1]);
        } else {
            if (tL1<tR1)
                dfs_rayHit(u->son[1]);
            if (tL0<tR0 && tL0<closest_hit.t_max)
                dfs_rayHit(u->son[0]);
        }
    }

    HitRecord rayHit(Ray _ray) {
        cnt = 0;
        ray = _ray;
        closest_hit = HitRecord(eps, INFINITY);
        dfs_rayHit(root);
        return closest_hit;
    }

    ~KDT() {
        delete buffer;
    }
};


#endif //_KDT_H
