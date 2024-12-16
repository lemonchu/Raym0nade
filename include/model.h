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
    vec3 shapeNormal, surfaceNormal, emission, baseColor, position;
    float roughness, metallic;
    HitInfo();
};

void getHitAllNormals(const Face& face, const vec3 &inDir, const vec3 &baryCoords,
                      vec3 &shapeNormal_to_compute, vec3 &surfaceNormal_to_compute_raw);
void getHitTexture(const Face& face, const vec3 &baryCoords, const vec3 &hit_dPdx, const vec3 &hit_dPdy, HitInfo &hitInfo_to_init);
void getHitInfo(const Face& face, const vec3 &inDir,
                HitInfo &hitInfo_to_init, const glm::vec3& hit_dPdx, const glm::vec3& hit_dPdy);

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
    [[nodiscard]] HitRecord rayHit(Ray ray) const;
    [[nodiscard]] bool rayHit_test(Ray ray, float aimDepth) const;
};

#endif // MODEL_H