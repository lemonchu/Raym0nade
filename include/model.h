#ifndef MODEL_H
#define MODEL_H

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "material.h"
#include "component.h"
#include "kdt.h"

struct HitInfo {
    vec3 shapeNormal, surfaceNormal, emission;
    vec4 diffuseColor;
    float t;
    HitInfo() : t(0), shapeNormal(vec3(0)), surfaceNormal(vec3(0)), emission(vec3(0)), diffuseColor(vec4(0)) {}
};

void getHitInfo(const Face& face, const glm::vec3& intersection, const vec3 &inDir, HitInfo &hitInfo);

class Model {
private:
    void processMesh(aiMesh *mesh, const glm::mat4 &nodeTransform);
    void processMaterial(std::string model_folder, const aiScene *scene);
    void processNode(aiNode *node, const aiScene *scene, const glm::mat4 &parentTransform);

public:
    std::vector<Material> materials;
    std::vector<Face> faces;
    std::vector<LightObject> lightObjects;
    KDT kdt;
    std::string model_path;

    Model();
    Model(const std::string &model_folder, const std::string &model_name);
    void checkEmissiveMaterials(const aiScene* scene);
    void rayHit(Ray ray, HitInfo &hitInfo) const;
    bool rayHit_test(Ray ray, float aimDepth) const;
    void rayCast(Ray ray, float &depth, const Face *&face) const;
};

#endif // MODEL_H