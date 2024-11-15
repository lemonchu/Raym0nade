#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <vector>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "geometry.h"
#include "texture.h"
#include "kdt.h"

class Model {
private:
    Assimp::Importer importer;
    const aiScene *scene;
public:
    std::vector<Texture> textures;
    std::vector<Triangle> triangles;
    KDT kdt;

    Model() {}

    int load(const char *file_name) {
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

        textures.resize(scene->mNumMaterials);
        for (unsigned int i = 0; i < scene->mNumMaterials; i++) {
            std::cout << "Loading material " << i << std::endl;
            aiMaterial *material = scene->mMaterials[i];

            for (unsigned int j = 0; j < AI_TEXTURE_TYPE_MAX; j++) {
        if (j!=1)continue;
                aiTextureType textureType = (aiTextureType) j;
                unsigned int numTextures = material->GetTextureCount(textureType);
                for (unsigned int k = 0; k < numTextures; k++) {
                    aiString path;
                    if (material->GetTexture(textureType, k, &path) == AI_SUCCESS) {
                        std::cout << "Texture path (" << textureType << "): " << path.C_Str() << std::endl;
                    }
                    std::string pathStr = (std::string)"fbx/" + path.C_Str();
#ifdef WIN32
                    std::replace(pathStr.begin(), pathStr.end(), '/', '\\');
#else
                    std::replace(pathStr.begin(), pathStr.end(), '\\', '/');
#endif
                    textures[i].loadImageFormFile(j,   pathStr);
                }
            }
        }
        std::cout << "Materials: " << textures.size() << std::endl;

        textures.reserve(scene->mNumMeshes);
        for (unsigned int i = 0; i < scene->mNumMeshes; i++) {
            printf("Loading mesh %d\n", i);
            aiMesh *mesh = scene->mMeshes[i];
            for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
                aiFace &face = mesh->mFaces[j];
                aiVector3D vertex[3] = {
                        mesh->mVertices[face.mIndices[0]],
                        mesh->mVertices[face.mIndices[1]],
                        mesh->mVertices[face.mIndices[2]]
                };
                aiVector3D texCoord[3] = {
                        mesh->mTextureCoords[0][face.mIndices[0]],
                        mesh->mTextureCoords[0][face.mIndices[1]],
                        mesh->mTextureCoords[0][face.mIndices[2]]
                };
                triangles.push_back({
                        {vec3(vertex[0].x, vertex[0].y, vertex[0].z),
                         vec3(vertex[1].x, vertex[1].y, vertex[1].z),
                         vec3(vertex[2].x, vertex[2].y, vertex[2].z)},
                        {glm::vec<2, float>(texCoord[0].x, texCoord[0].y),
                         glm::vec<2, float>(texCoord[1].x, texCoord[1].y),
                         glm::vec<2, float>(texCoord[2].x, texCoord[2].y)},
                        &textures[mesh->mMaterialIndex]
                });
            }
        }
        std::cout << "Faces: " << triangles.size() << std::endl;

        kdt.build(triangles);

        return 0;
    }

};

#endif // MODEL_H