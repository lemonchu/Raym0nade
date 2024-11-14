#ifndef MODEL_H
#define MODEL_H

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "geometry.h"
#include <iostream>
#include <vector>

class Model {
private:
    Assimp::Importer importer;
    const aiScene* scene;
public:
    std::vector<Triangle> triangles;
    int load(const char* file_name) {
        scene = importer.ReadFile(file_name,
                                  aiProcess_CalcTangentSpace |
                                  aiProcess_Triangulate |
                                  aiProcess_JoinIdenticalVertices |
                                  aiProcess_SortByPType |
                                  aiProcess_GenUVCoords);

        if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
            std::cerr << "Error loading model: " << importer.GetErrorString() << std::endl;
            return -1;
        }

        /*
        for (unsigned int i = 0; i < scene->mNumMaterials; i++) {
            aiMaterial* material = scene->mMaterials[i];

            for (unsigned int j = 0; j < AI_TEXTURE_TYPE_MAX; j++) {
                aiTextureType textureType = (aiTextureType)j;
                unsigned int numTextures = material->GetTextureCount(textureType);
                for (unsigned int k = 0; k < numTextures; k++) {
                    aiString path;
                    if (material->GetTexture(textureType, k, &path) == AI_SUCCESS) {
                        std::cout << "Texture path (" << textureType << "): " << path.C_Str() << std::endl;
                    }
                }
            }
        }
        */

        for (unsigned int i = 0; i < scene->mNumMeshes; i++) {
            aiMesh* mesh = scene->mMeshes[i];

            for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
                aiFace &face = mesh->mFaces[j];
                aiVector3D vertex[3] = {
                    mesh->mVertices[face.mIndices[0]],
                    mesh->mVertices[face.mIndices[1]],
                    mesh->mVertices[face.mIndices[2]]
                };
                triangles.push_back({{
                    glm::vec<3, float>(vertex[0].x, vertex[0].y, vertex[0].z),
                    glm::vec<3, float>(vertex[1].x, vertex[1].y, vertex[1].z),
                    glm::vec<3, float>(vertex[2].x, vertex[2].y, vertex[2].z)},
                    j
                });
            }
        }
        std::cout << "Faces: " << triangles.size() << std::endl;

        return 0;
    }

};

#endif // MODEL_H