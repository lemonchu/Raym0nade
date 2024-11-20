#ifndef MODEL_H
#define MODEL_H

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "geometry.h"
#include "material.h"
#include "kdt.h"



class Model {
private:
    const aiScene *scene;

    void processNode(aiNode *node, const aiScene *scene, const glm::mat4 &parentTransform);

public:
    std::vector<Material> materials;
    std::vector<Face> faces;
    std::vector<VertexData> vertexDatas;
    KDT kdt;
    std::string model_path;

    Model();
    Model(std::string model_folder, std::string model_name);
};

void checkEmissiveMaterials(const aiScene* scene);
void checkLightSources(const aiScene* scene);
glm::vec3 barycentric(const glm::vec3& A, const glm::vec3& B, const glm::vec3& C, const glm::vec3& P);

#endif // MODEL_H