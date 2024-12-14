#ifndef MODEL_H
#define MODEL_H

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "material.h"
#include "component.h"
#include "kdt.h"
#include <Python.h>

struct HitInfo {
    vec3 shapeNormal, surfaceNormal, emission;
    vec4 baseColor;
    float t, roughness, metallic;
    HitInfo() : t(0) {}
};

void getHitNormal(const Face& face, const vec3 &inDir, const vec3 &baryCoords, const vec2 &texUV,
                  vec3 &shapeNormal, vec3 &surfaceNormal);
void getHitInfo(const Face& face, const glm::vec3& intersection, const vec3 &inDir,
                HitInfo &hitInfo);

class Model {
private:
    void processMaterial(const std::string &model_folder, const aiScene *scene);
    void checkLightObject(Face *meshFaces, aiMesh *mesh, const Material &material);
    void processMesh(aiMesh *mesh, const glm::mat4 &nodeTransform);
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
    void rayHit(Ray ray, HitInfo &hitInfo) const;
    bool rayHit_test(Ray ray, float aimDepth) const;
    void rayCast(Ray ray, float &depth, const Face *&face) const;
};

#endif // MODEL_H