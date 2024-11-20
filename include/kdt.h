#ifndef KDT_H
#define KDT_H

#include "geometry.h"

struct KDT_Node {
    Box box;
    KDT_Node *son[2];
    Triangle *triL, *triR;
    KDT_Node();
};

class KDT {
    KDT_Node *buffer, *root;
    int cur;
public:
    KDT();

    KDT_Node *newNode();

    void dfs_build(KDT_Node *&u, Triangle *triL, Triangle *triR);

    void build(std::vector<Triangle> &triangles);

    HitRecord closest_hit;
    Ray ray;

    void dfs_rayHit(KDT_Node *u);

    HitRecord rayHit(Ray _ray);

    ~KDT();
};

#endif //_KDT_H
