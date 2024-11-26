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

LightFace::LightFace(vec3 position, vec3 normal, float power) :
        position(position), normal(normal), power(power) {}

LightObject::LightObject() : center(vec3(0)), color(vec3(0)), power(0) {}

unsigned int vertexCount(const aiScene *scene) {
    unsigned int cnt = 0;
    for (unsigned int i = 0; i < scene->mNumMeshes; i++)
        cnt += scene->mMeshes[i]->mNumVertices;
    return cnt;
}

void RandomNumberGenerator::Init(const std::vector<float>& distribution) {
    prefixSums.resize(distribution.size());
    prefixSums[0] = distribution[0];
    for (size_t i = 1; i < distribution.size(); ++i) {
        prefixSums[i] = prefixSums[i - 1] + distribution[i];
    }
}

int RandomNumberGenerator::operator()(std::mt19937 &gen) const {
    std::uniform_real_distribution<float> dist(0.0f, prefixSums.back());
    float randomValue = dist(gen);
    return std::lower_bound(prefixSums.begin(), prefixSums.end(), randomValue) - prefixSums.begin();
}

const float powerAlpha = 0.25f;

void Model::processMesh(aiMesh *mesh, const glm::mat4 &nodeTransform) {

    unsigned int offset = vertexDatas.size();
    for (unsigned int j = 0; j < mesh->mNumVertices; j++) {
        aiVector3D normal = mesh->mNormals[j];
        aiVector3D texCoord = mesh->mTextureCoords[0][j];
        vertexDatas.emplace_back(
                vec2(texCoord.x, texCoord.y),
                vec3(0)
        );
    }

    auto &material = materials[mesh->mMaterialIndex];
    VertexData *vData = &vertexDatas[offset];
    offset = faces.size();
    for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
        aiFace &face0 = mesh->mFaces[j];

        // Apply transformation to vertices
        glm::vec4 transformedVertex[3];
        for (int k = 0; k < 3; ++k) {
            aiVector3D vertex = mesh->mVertices[face0.mIndices[k]];
            transformedVertex[k] = nodeTransform * glm::vec4(vertex.x, vertex.y, vertex.z, 1.0f);
        }

        faces.push_back({
                {vec3(transformedVertex[0]),
                 vec3(transformedVertex[1]),
                 vec3(transformedVertex[2])},
                {&vData[face0.mIndices[0]],
                 &vData[face0.mIndices[1]],
                 &vData[face0.mIndices[2]]},
                &material
        });

        const Face &face = faces.back();

        vec3 normal = glm::normalize(glm::cross(face.v[1] - face.v[0], face.v[2] - face.v[0]));
        vData[face0.mIndices[0]].normal += normal;
        vData[face0.mIndices[1]].normal += normal;
        vData[face0.mIndices[2]].normal += normal;
    }

    for (unsigned int j = 0; j < mesh->mNumVertices; j++)
        vData[j].normal = glm::normalize(vData[j].normal);

    Face *meshFaces = &faces[offset];
    if (material.isEmissionEnabled && !material.texture[TextureIdForEmission].empty()) {

        lightObjects.emplace_back();
        auto &lightObject = lightObjects.back();

        vec3 emission = material.emission;
        lightObject.color = glm::normalize(emission);

        float totalColorPower = 0.0f;

        for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
            auto &face = meshFaces[j];

            vec3 centroid = (face.v[0] + face.v[1] + face.v[2]) / 3.0f;
            vec2 uv = (face.data[0]->uv + face.data[1]->uv + face.data[2]->uv) / 3.0f;
            vec4 textureColor = material.getImage(TextureIdForEmission).get(uv[0], uv[1]);
            float colorPower = glm::length(vec3(textureColor[0], textureColor[1], textureColor[2]));

            if (colorPower == 0.0f)
                continue;

            float area = glm::length(glm::cross(face.v[1] - face.v[0], face.v[2] - face.v[0])) / 2.0f;
            vec3 normal = normalize(face.data[0]->normal + face.data[1]->normal + face.data[2]->normal);
            lightObject.lightFaces.emplace_back(centroid, normal, area * colorPower);
            totalColorPower += area * colorPower;
        }

        float powerDensity = glm::length(emission) * pow(totalColorPower, powerAlpha-1);
        std::vector<float> faceWeights;
        faceWeights.resize(lightObject.lightFaces.size());
        for (unsigned int j = 0; j < lightObject.lightFaces.size(); j++) {
            LightFace &lightFace = lightObject.lightFaces[j];
            lightFace.power *= powerDensity;
            faceWeights[j] = lightFace.power;
            lightObject.power += lightFace.power;
            lightObject.center += lightFace.position * lightFace.power;
        }

        if (lightObject.power == 0.0f)
            lightObjects.pop_back();
        else {
                lightObject.faceDist.Init(faceWeights);
                lightObject.center /= lightObject.power;

                std::cout << "Light object with power: " << lightObject.power
                          << ", color:" << lightObject.color.x << ", " << lightObject.color.y << ", " << lightObject.color.z
                          << ", position:" << lightObject.center.x << ", " << lightObject.center.y << ", "
                          << lightObject.center.z << std::endl;
        }
    }
}

void Model::processNode(aiNode *node, const aiScene *scene, const glm::mat4 &parentTransform) {
    glm::mat4 nodeTransform = parentTransform * aiMatrix4x4ToGlm(node->mTransformation);

    for (unsigned int i = 0; i < node->mNumMeshes; i++) {
        aiMesh *mesh = scene->mMeshes[node->mMeshes[i]];
        processMesh(mesh, nodeTransform);
    }

    for (unsigned int i = 0; i < node->mNumChildren; i++) {
        processNode(node->mChildren[i], scene, nodeTransform);
    }
}

void Model::processMaterial(std::string model_folder, const aiScene *scene) {

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
                    std::cout << "- Texture path (" << textureType << "): " << urlDecode(path.C_Str()) << std::endl;
                }
                std::string pathStr = (std::string) model_folder + path.C_Str();
#ifdef WIN32
                std::replace(pathStr.begin(), pathStr.end(), '/', '\\');
#else
                std::replace(pathStr.begin(), pathStr.end(), '\\', '/');
#endif
                pathStr = urlDecode(pathStr);
                materials[i].loadImageFromFile(j, pathStr);
            }
        }
        materials[i].loadMaterialProperties(material);
    }
    std::cout << "Materials: " << materials.size() << std::endl;
}

Model::Model() = default;

Model::Model(const std::string &model_folder, const std::string &model_name) {

    Assimp::Importer importer;
    model_path = model_folder + model_name;
    const aiScene* scene = importer.ReadFile(model_path.c_str(),
                                             aiProcessPreset_TargetRealtime_Quality);

    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "Error loading model: " << importer.GetErrorString() << std::endl;
        return ;
    }
    unsigned int vertexCnt = vertexCount(scene);
    std::cout << "Vertices: " << vertexCnt << std::endl;
    vertexDatas.reserve(vertexCnt);

    processMaterial(model_folder, scene);

    const glm::mat4 identity = glm::mat4(1.0f);
    processNode(scene->mRootNode, scene, identity);

    checkEmissiveMaterials(scene);
    checkLightSources(scene);

    importer.FreeScene();

    std::cout << "Faces: " << faces.size() << std::endl;

    kdt.build(faces);
}

void Model::checkEmissiveMaterials(const aiScene* scene) {
    for (unsigned int i = 0; i < scene->mNumMaterials; i++) {
        aiMaterial* material = scene->mMaterials[i];
        aiColor4D emissiveColor(0.0f,0.0f,0.0f,0.0f);
        if (AI_SUCCESS == material->Get(AI_MATKEY_COLOR_EMISSIVE, emissiveColor)) {
            materials[i].emission = glm::vec3(emissiveColor.r, emissiveColor.g, emissiveColor.b);
            if (emissiveColor.r > 0.0f || emissiveColor.g > 0.0f || emissiveColor.b > 0.0f) {
                std::cout << "Material " << i << " has emissive color: "
                          << emissiveColor.r << ", " << emissiveColor.g << ", " << emissiveColor.b << ", " << emissiveColor.a << std::endl;
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