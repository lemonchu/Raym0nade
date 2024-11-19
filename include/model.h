#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <vector>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "geometry.h"
#include "texture.h"
#include "kdt.h"

glm::mat4 aiMatrix4x4ToGlm(const aiMatrix4x4 &from) {
    glm::mat4 to;
    to[0][0] = from.a1;
    to[0][1] = from.b1;
    to[0][2] = from.c1;
    to[0][3] = from.d1;
    to[1][0] = from.a2;
    to[1][1] = from.b2;
    to[1][2] = from.c2;
    to[1][3] = from.d2;
    to[2][0] = from.a3;
    to[2][1] = from.b3;
    to[2][2] = from.c3;
    to[2][3] = from.d3;
    to[3][0] = from.a4;
    to[3][1] = from.b4;
    to[3][2] = from.c4;
    to[3][3] = from.d4;
    return to;
}

class Model {
private:
    Assimp::Importer importer;
    const aiScene *scene;

    void processNode(aiNode *node, const aiScene *scene, const glm::mat4 &parentTransform) {
        glm::mat4 nodeTransform = parentTransform * aiMatrix4x4ToGlm(node->mTransformation);

        for (unsigned int i = 0; i < node->mNumMeshes; i++) {
            aiMesh *mesh = scene->mMeshes[node->mMeshes[i]];

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

                // Apply transformation to vertices
                glm::vec4 transformedVertex[3];
                for (int k = 0; k < 3; ++k) {
                    transformedVertex[k] = nodeTransform * glm::vec4(vertex[k].x, vertex[k].y, vertex[k].z, 1.0f);
                }

                triangles.push_back({{vec3(transformedVertex[0]), vec3(transformedVertex[1]),
                                      vec3(transformedVertex[2])},
                                     {glm::vec<2, float>(texCoord[0].x, texCoord[0].y),
                                      glm::vec<2, float>(texCoord[1].x, texCoord[1].y),
                                      glm::vec<2, float>(texCoord[2].x, texCoord[2].y)},
                                     &textures[mesh->mMaterialIndex]
                                    });
            }
        }

        for (unsigned int i = 0; i < node->mNumChildren; i++) {
            processNode(node->mChildren[i], scene, nodeTransform);
        }
    }

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
                aiTextureType textureType = (aiTextureType) j;
                unsigned int numTextures = material->GetTextureCount(textureType);
                for (unsigned int k = 0; k < numTextures; k++) {
                    aiString path;
                    if (material->GetTexture(textureType, k, &path) == AI_SUCCESS) {
                        std::cout << "Texture path (" << textureType << "): " << path.C_Str() << std::endl;
                    }
                    std::string pathStr = (std::string) "fbx/" + path.C_Str();
#ifdef WIN32
                    std::replace(pathStr.begin(), pathStr.end(), '/', '\\');
#else
                    std::replace(pathStr.begin(), pathStr.end(), '\\', '/');
#endif
                    textures[i].loadImageFromFile(j, pathStr);
                }
            }
        }
        std::cout << "Materials: " << textures.size() << std::endl;

        glm::mat4 identity = glm::mat4(1.0f);
        processNode(scene->mRootNode, scene, identity);

        std::cout << "Faces: " << triangles.size() << std::endl;

        kdt.build(triangles);

        return 0;
    }

};

#endif // MODEL_H