#ifndef KDT_H
#define KDT_H

#include "component.h"
#include "geometry.h"

struct BVH_Node {
    Box box;
    int faceL, faceR;
    BVH_Node();
};

class BVH {
    BVH_Node *node;
    Face *faces;
public:

    BVH();

    void dfs_build(int u, int faceL, int faceR);

    void build(std::vector<Face> &faces);

    void dfs_rayHit(int u, const Ray &ray, HitRecord &closest_hit) const;

    void rayHit(const Ray &ray, HitRecord &closest_hit) const;

    ~BVH();
};

#endif //KDT_H