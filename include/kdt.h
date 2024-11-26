#ifndef KDT_H
#define KDT_H

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

    [[nodiscard]] HitRecord rayHit(const Ray &ray, float t_max = INFINITY) const;

    ~KDT();
};

#endif //_KDT_H
