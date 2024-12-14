#ifndef KDT_H
#define KDT_H

#include "component.h"
#include "geometry.h"

struct KDT_Node {
    Box box;
    int faceL, faceR;
    KDT_Node();
};

class KDT {
    KDT_Node *node;
    Face *faces;
public:

    KDT();

    void dfs_build(int u, int faceL, int faceR);

    void build(std::vector<Face> &faces);

    void dfs_rayHit(int u, const Ray &ray, HitRecord &closest_hit) const;

    void rayHit(const Ray &ray, HitRecord &closest_hit) const;

    ~KDT();
};

#endif //KDT_H
