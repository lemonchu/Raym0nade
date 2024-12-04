#ifndef KDT_H
#define KDT_H

#include "component.h"
#include "geometry.h"

struct KDT_Node {
    Box box;
    KDT_Node *son[2];
    Face *faceL, *faceR;
    KDT_Node();
};

class KDT {
    KDT_Node *buffer, *root;
    int cur;
public:
    KDT();

    KDT_Node *newNode();

    void dfs_build(KDT_Node *&u, Face *faceL, Face *faceR);

    void build(std::vector<Face> &faces);

    void dfs_rayHit(KDT_Node *u, const Ray &ray, HitRecord &closest_hit) const;

    void rayHit(const Ray &ray, HitRecord &closest_hit) const;

    ~KDT();
};

#endif //KDT_H
