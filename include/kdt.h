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

    HitRecord closest_hit;
    Ray ray;

    void dfs_rayHit(KDT_Node *u);

    HitRecord rayHit(Ray _ray);

    ~KDT();
};

#endif //_KDT_H
