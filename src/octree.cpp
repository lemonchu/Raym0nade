#include <iostream>
#include "octree.h"

OctreeNode::OctreeNode() {
    child[0][0][0] = child[0][0][1] = child[0][1][0] = child[0][1][1] =
    child[1][0][0] = child[1][0][1] = child[1][1][0] = child[1][1][1] = nullptr;
    probe = nullptr;
    weight = 0.0f;
}

const float ClusterRatio = 0.8f;

void initProbe(OctreeNode *u, Importron* L, Importron* R, float totalWeight) {
    Probe &probe = *(u->probe = new Probe);
    vec3 normal = vec3(0.0f);
    for (Importron *ptr = L; ptr != R; ptr++)
        normal += ptr->normal * ptr->weight;
    normal = normalize(normal);
    float Edot = 0.0f;
    for (Importron *ptr = L; ptr != R; ptr++)
        Edot += std::max(glm::dot(normal, ptr->normal), 0.0f) * ptr->weight;
    Edot /= totalWeight;
    for (Importron *ptr = L; ptr != R; ptr++)
        if (glm::dot(normal, ptr->normal) > Edot * ClusterRatio)
            probe.normal += ptr->normal * ptr->weight;
    probe.normal = normalize(probe.normal);
    getTangentSpace(probe.normal, probe.tangent, probe.bitangent);
}

Importron* sort01(Importron* L, Importron* R, int axis, float M) {
    R--;
    while(true) {
        while(L <= R && L->pos[axis]<=M)
            L++;
        while(L <= R && R->pos[axis]>=M)
            R--;
        if (L > R)
            break;
        std::swap(*L, *R);
        L++;
        R--;
    }
    return L;
}

Octree::Octree() {
    root = nullptr;
}

void Octree::dfs_build(OctreeNode *&u, Importron* L, Importron* R, vec3 v0, vec3 v1, int depth) {
    float totalWeight = 0.0f;
    for (Importron* ptr = L; ptr != R; ptr++)
        totalWeight += ptr->weight;

    if (totalWeight < minWeight)
        return ;

    u = new OctreeNode();
    u->weight = totalWeight;
    if (depth == maxDepth) {
        initProbe(u, L, R, totalWeight);
        return ;
    }
    vec3 vM = (v0 + v1) / 2.0f;
    u->vM = vM;

    Importron* P[9];
    P[0] = L;
    P[8] = R;
    P[4] = sort01(P[0], P[8], 0, vM[0]);
    P[2] = sort01(P[0], P[4], 1, vM[1]);
    P[6] = sort01(P[4], P[8], 1, vM[1]);
    P[1] = sort01(P[0], P[2], 2, vM[2]);
    P[3] = sort01(P[2], P[4], 2, vM[2]);
    P[5] = sort01(P[4], P[6], 2, vM[2]);
    P[7] = sort01(P[6], P[8], 2, vM[2]);

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                vec3 v0_ = vec3(i, j, k) * (v1 - v0) / 2.0f + v0;
                vec3 v1_ = vec3(i + 1, j + 1, k + 1) * (v1 - v0) / 2.0f + v0;
                dfs_build(u->child[i][j][k], P[i*4+j*2+k], P[i*4+j*2+k+1], v0_, v1_, depth + 1);
                if (u->child[i][j][k] != nullptr)
                    u->weight -= u->child[i][j][k]->weight;
            }

    if (u->weight > minWeight && depth >= minDepth)
        initProbe(u, L, R, totalWeight);
    else u->weight = 0;
}

void Octree::build(std::vector<Importron> &importrons) {
    dfs_build(root, &importrons.front(), &importrons.back()+1, vec3(-coordLimit), vec3(coordLimit), 0);
}

Probe* Octree::dfs_query(OctreeNode *u, const vec3 &pos) const{
    if (u == nullptr)
        return nullptr;
    Probe *ret = dfs_query(u->child[pos[0]>u->vM[0]][pos[1]>u->vM[1]][pos[2]>u->vM[2]], pos);
    return (ret == nullptr) ? u->probe : ret;
}

Probe* Octree::query(const vec3 &pos) const{
    return dfs_query(root, pos);
}

void Octree::dfs_delete(OctreeNode *&u) {
    if (u == nullptr)
        return ;
    delete u->probe;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                dfs_delete(u->child[i][j][k]);
    delete u;
}

Octree::~Octree() {
    dfs_delete(root);
}