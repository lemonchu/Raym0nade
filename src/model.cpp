#include <iostream>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "model.h"

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

unsigned int vertexCount(const aiScene *scene) {
    unsigned int cnt = 0;
    for (unsigned int i = 0; i < scene->mNumMeshes; i++)
        cnt += scene->mMeshes[i]->mNumVertices;
    return cnt;
}

void Model::processNode(aiNode *node, const aiScene *scene, const glm::mat4 &parentTransform) {
    glm::mat4 nodeTransform = parentTransform * aiMatrix4x4ToGlm(node->mTransformation);

    for (unsigned int i = 0; i < node->mNumMeshes; i++) {
        aiMesh *mesh = scene->mMeshes[node->mMeshes[i]];
        unsigned int offset = vertexDatas.size();
        bool hasVertexColor = mesh->HasVertexColors(0);
        for (unsigned int j = 0; j < mesh->mNumVertices; j++) {
            aiVector3D normal = mesh->mNormals[j];
            aiVector3D texCoord = mesh->mTextureCoords[0][j];
            vertexDatas.emplace_back(
                    vec2(texCoord.x, texCoord.y),
                    vec3(normal.x, normal.y, normal.z)
            );
        }
        VertexData *vData = &vertexDatas[offset];
        for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
            aiFace &face = mesh->mFaces[j];

            // Apply transformation to vertices
            glm::vec4 transformedVertex[3];
            for (int k = 0; k < 3; ++k) {
                aiVector3D vertex = mesh->mVertices[face.mIndices[k]];
                transformedVertex[k] = nodeTransform * glm::vec4(vertex.x, vertex.y, vertex.z, 1.0f);
            }

            faces.push_back({
                {vec3(transformedVertex[0]),
                 vec3(transformedVertex[1]),
                 vec3(transformedVertex[2])},
                 {&vData[face.mIndices[0]],
                 &vData[face.mIndices[1]],
                 &vData[face.mIndices[2]]},
                 &materials[mesh->mMaterialIndex]
            });
        }
    }

    for (unsigned int i = 0; i < node->mNumChildren; i++) {
        processNode(node->mChildren[i], scene, nodeTransform);
    }
}

Model::Model() = default;

Model::Model(std::string model_folder, std::string model_name) {

    Assimp::Importer importer;
    model_path = model_folder + model_name;
    scene = importer.ReadFile(model_path.c_str(),
                              aiProcess_CalcTangentSpace |
                              aiProcess_Triangulate |
                              aiProcess_JoinIdenticalVertices |
                              aiProcess_SortByPType |
                              aiProcess_GenUVCoords |
                              aiProcess_GenNormals);

    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "Error loading model: " << importer.GetErrorString() << std::endl;
        return ;
    }
    unsigned int vertexCnt = vertexCount(scene);
    std::cout << "Vertices: " << vertexCnt << std::endl;
    vertexDatas.reserve(vertexCnt);

    materials.resize(scene->mNumMaterials);
    for (unsigned int i = 0; i < scene->mNumMaterials; i++) {
        std::cout << "Loading material " << i << std::endl;
        aiMaterial *material = scene->mMaterials[i];

        for (unsigned int j = 0; j < AI_TEXTURE_TYPE_MAX; j++) {
            aiTextureType textureType = (aiTextureType) j;
            unsigned int numTextures = material->GetTextureCount(textureType);
            for (unsigned int k = 0; k < numTextures; k++) {
                aiString path;
                if (material->GetTexture(textureType, k, &path) == AI_SUCCESS) {
                    std::cout << "Material path (" << textureType << "): " << path.C_Str() << std::endl;
                }
                std::string pathStr = (std::string) model_folder + path.C_Str();
#ifdef WIN32
                std::replace(pathStr.begin(), pathStr.end(), '/', '\\');
#else
                std::replace(pathStr.begin(), pathStr.end(), '\\', '/');
#endif
                materials[i].loadImageFromFile(j, pathStr);
            }
        }
        materials[i].loadMaterialProperties(material);
    }
    std::cout << "Materials: " << materials.size() << std::endl;

    const glm::mat4 identity = glm::mat4(1.0f);
    processNode(scene->mRootNode, scene, identity);

    checkEmissiveMaterials(scene);
    checkLightSources(scene);

    importer.FreeScene();

    std::cout << "Faces: " << faces.size() << std::endl;

    kdt.build(faces);
}

void checkEmissiveMaterials(const aiScene* scene) {
    for (unsigned int i = 0; i < scene->mNumMaterials; i++) {
        aiMaterial* material = scene->mMaterials[i];
        aiColor3D emissiveColor(0.0f, 0.0f, 0.0f);
        if (AI_SUCCESS == material->Get(AI_MATKEY_COLOR_EMISSIVE, emissiveColor)) {
            if (emissiveColor.r > 0.0f || emissiveColor.g > 0.0f || emissiveColor.b > 0.0f) {
                std::cout << "Material " << i << " has emissive color: "
                          << emissiveColor.r << ", " << emissiveColor.g << ", " << emissiveColor.b << std::endl;
            }
        }
    }
}

void checkLightSources(const aiScene* scene) {
    for (unsigned int i = 0; i < scene->mNumLights; i++) {
        aiLight* light = scene->mLights[i];
        std::cout << "Light " << i << " of type " << light->mType << " with color: "
                  << light->mColorDiffuse.r << ", " << light->mColorDiffuse.g << ", " << light->mColorDiffuse.b << std::endl;
    }
}