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

#endif // MODEL_H