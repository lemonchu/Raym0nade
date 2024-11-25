#ifndef MODEL_H
#define MODEL_H

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <random>
#include "geometry.h"
#include "material.h"
#include "kdt.h"

struct LightFace {
    vec3 position, normal;
    float power;
    LightFace(vec3 position, vec3 normal, float power);
};

class LightObject {
public:
    vec3 center, color;
    float power;
    unsigned long long touch, total;
    std::vector<LightFace> lightFaces;
    std::discrete_distribution<int> faceDist;
    LightObject();
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
    Model(std::string model_folder, std::string model_name);
    void checkEmissiveMaterials(const aiScene* scene);
};

void checkEmissiveMaterials(const aiScene* scene);
void checkLightSources(const aiScene* scene);
glm::vec3 barycentric(const glm::vec3& A, const glm::vec3& B, const glm::vec3& C, const glm::vec3& P);

#endif // MODEL_H