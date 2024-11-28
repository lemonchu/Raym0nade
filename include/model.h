#ifndef MODEL_H
#define MODEL_H

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "material.h"
#include "component.h"
#include "kdt.h"

struct HitInfo {
    float t;
    vec3 shapeNormal, surfaceNormal, emission;
    vec4 diffuseColor;
    HitInfo() : t(0), shapeNormal(vec3(0)), surfaceNormal(vec3(0)), emission(vec3(0)), diffuseColor(vec4(0)) {}
};

class Model {
private:
    void processMesh(aiMesh *mesh, const glm::mat4 &nodeTransform);
    void processMaterial(std::string model_folder, const aiScene *scene);
    void processNode(aiNode *node, const aiScene *scene, const glm::mat4 &parentTransform);

public:
    std::vector<Material> materials;
    std::vector<Face> faces;
    std::vector<VertexData> vertexDatas;
    std::vector<LightObject> lightObjects;
    KDT kdt;
    std::string model_path;

    Model();
    Model(const std::string &model_folder, const std::string &model_name);
    void checkEmissiveMaterials(const aiScene* scene);
    HitInfo rayHit(Ray ray) const;
    bool rayHit_test(Ray ray, float aimDepth) const;
};

#endif // MODEL_H