#ifndef OCTREE_H
#define OCTREE_H

#include <random>
#include "geometry.h"
#include "sampling.h"

struct Importron {
    vec3 pos;
    union {
        vec3 normal;
        vec3 direction;
    };
    union {
        float weight;
        float power;
    };
};

struct OctreeNode {
    OctreeNode *child[2][2][2];
    vec3 vM;
    Probe *probe;
    float weight;
    OctreeNode();
};

class Octree {
private:
    void dfs_build(OctreeNode *&u, Importron* L, Importron* R, vec3 v0, vec3 v1, int depth = 0);
    Probe* dfs_query(OctreeNode *u, const vec3 &pos) const;
    void dfs_delete(OctreeNode *&u);
public:
    OctreeNode *root;
    float coordLimit, minWeight;
    unsigned int minDepth, maxDepth;
    Octree();
    void build(std::vector<Importron> &importrons);
    [[nodiscard]] Probe* query(const vec3 &pos) const;
    ~Octree();
};

#endif // OCTREE_H
