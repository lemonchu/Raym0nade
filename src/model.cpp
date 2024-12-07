#include <iostream>
#include "model.h"

unsigned int vertexCount(const aiScene *scene) {
    unsigned int cnt = 0;
    for (int i = 0; i < scene->mNumMeshes; i++)
        cnt += scene->mMeshes[i]->mNumVertices;
    return cnt;
}

glm::mat4 aiMatrix4x4ToGlm(const aiMatrix4x4 &from) {
    glm::mat4 to;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            to[i][j] = from[j][i];
    return to;
}

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
                                &material,
                                nullptr
                        });

        const Face &face = faces.back();

        vec3 normal = normalize(glm::cross(face.v[1] - face.v[0], face.v[2] - face.v[0]));
        vData[face0.mIndices[0]].normal += normal;
        vData[face0.mIndices[1]].normal += normal;
        vData[face0.mIndices[2]].normal += normal;
    }

    for (unsigned int j = 0; j < mesh->mNumVertices; j++)
        vData[j].normal = normalize(vData[j].normal);

    Face *meshFaces = &faces[offset];
    if (material.isEmissionEnabled && !material.texture[aiTextureType_EMISSIVE].empty()) {

        lightObjects.emplace_back();
        auto &lightObject = lightObjects.back();

        float totalPower = 0.0f;
        vec3 color = vec3(0.0f);

        for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
            auto &face = meshFaces[j];

            vec3 centroid = (face.v[0] + face.v[1] + face.v[2]) / 3.0f;
            vec2 uv = (face.data[0]->uv + face.data[1]->uv + face.data[2]->uv) / 3.0f;
            vec3 textureColor = material.getEmissiveColor(uv[0], uv[1]);

            float colorPower = glm::length(textureColor);
            if (colorPower == 0.0f)
                continue;
            vec3 normal = cross(face.v[1] - face.v[0], face.v[2] - face.v[0]);
            float area = length(normal) / 2.0f;
            if (area < eps_zero)
                continue;
            normal = normalize(normal);
            color += textureColor * area;
            totalPower += area * colorPower;
            lightObject.lightFaces.emplace_back(centroid, normal, area * colorPower);
            face.lightObject = &lightObject;
        }
        if (lightObject.lightFaces.empty()) {
            lightObjects.pop_back();
        } else {
            static const float powerAlpha = 0.75f;
            lightObject.powerDensity = pow(totalPower, powerAlpha-1);
            lightObject.color = normalize(color);
            std::vector<float> faceWeights;
            faceWeights.resize(lightObject.lightFaces.size());
            for (unsigned int j = 0; j < lightObject.lightFaces.size(); j++) {
                LightFace &lightFace = lightObject.lightFaces[j];
                lightFace.power *= lightObject.powerDensity;
                faceWeights[j] = lightFace.power;
                lightObject.power += lightFace.power;
                lightObject.center += lightFace.position * lightFace.power;
            }

            lightObject.faceDist.Init(faceWeights);
            lightObject.center /= lightObject.power;

            std::cout << "Light object with power: " << lightObject.power
                      << ", power density: " << lightObject.powerDensity
                      << ", color:" << lightObject.color.x << ", " << lightObject.color.y << ", " << lightObject.color.z
                      << ", position:" << lightObject.center.x << ", " << lightObject.center.y << ", "
                      << lightObject.center.z << std::endl;
        }
    }
}

void Model::processNode(aiNode *node, const aiScene *scene, const glm::mat4 &parentTransform) {
    glm::mat4 nodeTransform = parentTransform * aiMatrix4x4ToGlm(node->mTransformation);

    for (int i = 0; i < node->mNumMeshes; i++) {
        aiMesh *mesh = scene->mMeshes[node->mMeshes[i]];
        processMesh(mesh, nodeTransform);
    }

    for (int i = 0; i < node->mNumChildren; i++) {
        processNode(node->mChildren[i], scene, nodeTransform);
    }
}

void Model::processMaterial(const std::string &model_folder, const aiScene *scene) {

    materials.resize(scene->mNumMaterials);
    for (int i = 0; i < scene->mNumMaterials; i++) {
        std::cout << "Loading material " << i << std::endl;
        aiMaterial *material = scene->mMaterials[i];

        for (int j = 0; j < AI_TEXTURE_TYPE_MAX; j++) {
            aiTextureType textureType = (aiTextureType) j;
            unsigned int numTextures = material->GetTextureCount(textureType);
            for (int k = 0; k < numTextures; k++) {
                aiString path;
                if (material->GetTexture(textureType, k, &path) == AI_SUCCESS) {
                    std::cout << "- Texture path (" << textureType << "): " << urlDecode(path.C_Str()) << std::endl;
                }
                std::string pathStr = model_folder + path.C_Str();
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

    // Finalize Python
    if (Py_IsInitialized()) {
        Py_Finalize();
    }
}

Model::Model() = default;

Model::Model(const std::string &model_folder, const std::string &model_name) {

    Assimp::Importer importer;
    model_path = model_folder + model_name;
    const aiScene* scene = importer.ReadFile(model_path.c_str(),
                                             aiProcess_Triangulate |
                                             aiProcess_JoinIdenticalVertices);

    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "Error loading model: " << importer.GetErrorString() << std::endl;
        return ;
    }
    unsigned int vertexCnt = vertexCount(scene);
    std::cout << "Vertices: " << vertexCnt << std::endl;
    vertexDatas.reserve(vertexCnt);
    lightObjects.reserve(scene->mNumMeshes);

    processMaterial(model_folder, scene);

    const glm::mat4 identity = glm::mat4(1.0f);
    processNode(scene->mRootNode, scene, identity);

    checkEmissiveMaterials(scene);
    //checkLightSources(scene);

    importer.FreeScene();

    std::cout << "Faces: " << faces.size() << std::endl;

    kdt.build(faces);
}

void Model::checkEmissiveMaterials(const aiScene* scene) {
    for (int i = 0; i < scene->mNumMaterials; i++) {
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
    for (int i = 0; i < scene->mNumLights; i++) {
        aiLight* light = scene->mLights[i];
        std::cout << "Light " << i << " of type " << light->mType << " with color: "
                  << light->mColorDiffuse.r << ", " << light->mColorDiffuse.g << ", " << light->mColorDiffuse.b << std::endl;
    }
}

bool TransparentTest(const Ray &ray, const HitRecord &hit) {
    const Face& face = *hit.face;
    const Material& material = *face.material;
    if (!material.hasTransparentPart)
        return false;
    const vec3 intersection = ray.origin + ray.direction * hit.t_max;
    const vec3 baryCoords = barycentric(face.v[0], face.v[1], face.v[2], intersection);
    vec2 texUV =
            baryCoords[0] * face.data[0]->uv
            + baryCoords[1] * face.data[1]->uv
            + baryCoords[2] * face.data[2]->uv;
    vec4 diffuseColor = material.getDiffuseColor(texUV[0], texUV[1]);
    return diffuseColor[3] == 0.0f;
}

void reverseCheck(vec3 &v, const vec3 &Dir) {
    if (dot(v, Dir) < 0.0f)
        v *= -1.0f;
}

void getHitInfo(const Face& face, const vec3& intersection, const vec3 &inDir, HitInfo &hitInfo) {
    const vec3 baryCoords = barycentric(face.v[0], face.v[1], face.v[2], intersection);
    const Material& material = *face.material;
    vec2 texUV =
            baryCoords[0] * face.data[0]->uv
            + baryCoords[1] * face.data[1]->uv
            + baryCoords[2] * face.data[2]->uv;

    hitInfo.diffuseColor = material.getDiffuseColor(texUV[0], texUV[1]);

    vec3 edge1 = face.v[1] - face.v[0];
    vec3 edge2 = face.v[2] - face.v[0];
    hitInfo.shapeNormal = normalize(cross(edge1, edge2));
    reverseCheck(hitInfo.shapeNormal, -inDir);

    vec3 normalV0 = face.data[0]->normal;
    vec3 normalV1 = face.data[1]->normal;
    vec3 normalV2 = face.data[2]->normal;
    reverseCheck(normalV0, hitInfo.shapeNormal);
    reverseCheck(normalV1, hitInfo.shapeNormal);
    reverseCheck(normalV2, hitInfo.shapeNormal);
    static const float threshold = 0.8f;
    hitInfo.surfaceNormal = normalize(
            baryCoords[0] * (dot(normalV0, hitInfo.shapeNormal) > threshold ? normalV0 : hitInfo.shapeNormal)
            + baryCoords[1] * (dot(normalV1, hitInfo.shapeNormal) > threshold ? normalV1 : hitInfo.shapeNormal)
            + baryCoords[2] * (dot(normalV2, hitInfo.shapeNormal) > threshold ? normalV2 : hitInfo.shapeNormal)
    );

    vec2 deltaUV1 = face.data[1]->uv - face.data[0]->uv;
    vec2 deltaUV2 = face.data[2]->uv - face.data[0]->uv;
    float f = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);
    vec3 tangent = f * (deltaUV2.y * edge1 - deltaUV1.y * edge2);
    vec3 bitangent = f * (-deltaUV2.x * edge1 + deltaUV1.x * edge2);
    tangent = normalize(tangent - hitInfo.shapeNormal * dot(hitInfo.shapeNormal, tangent));
    bitangent = normalize(bitangent - hitInfo.shapeNormal * dot(hitInfo.shapeNormal, bitangent) - tangent * dot(tangent, bitangent));
    vec3 normalMap = material.getNormal(texUV[0], texUV[1]);
    normalMap.z = 1.0f;
    hitInfo.surfaceNormal = normalize(
            tangent * normalMap.x
            + bitangent * normalMap.y
            + hitInfo.surfaceNormal * normalMap.z
    );
    if (isnan(hitInfo.surfaceNormal))
        hitInfo.surfaceNormal = hitInfo.shapeNormal;
    hitInfo.emission = (face.lightObject) ?
            material.getEmissiveColor(texUV[0], texUV[1]) * face.lightObject->powerDensity
            : vec3(0.0f);
}

const int maxRayDepth_hit = 8;

bool Model::rayHit_test(Ray ray, float aimDepth) const {
    HitRecord hit(eps_zero, aimDepth + eps_zero);
    for (int T = 0; T < maxRayDepth_hit; T++) {
        kdt.rayHit(ray, hit);
        if (hit.t_max > aimDepth)
            return false;
        if (!TransparentTest(ray, hit))
            return true;
        hit = HitRecord(hit.t_max + eps_zero, aimDepth + eps_zero);
    }
    return true;
}

void Model::rayHit(Ray ray, HitInfo &hitInfo) const {
    HitRecord hit;
    for (int T = 0; T < maxRayDepth_hit; T++) {
        kdt.rayHit(ray, hit);
        if (hit.t_max == INFINITY) {
            hitInfo.t = INFINITY;
            break;
        }
        vec3 intersection = ray.origin + ray.direction * hit.t_max;
        getHitInfo(*hit.face, intersection, ray.direction, hitInfo);
        hitInfo.t = hit.t_max;
        if (hitInfo.diffuseColor[3] > 0.0f)
            return;
        hit = HitRecord(hit.t_max + eps_zero, INFINITY);
    }
}

void Model::rayCast(Ray ray, float &depth, const Face *&face) const {
    HitRecord hit(eps_zero, INFINITY);
    for (int T = 0; T < maxRayDepth_hit; T++) {
        kdt.rayHit(ray, hit);
        if (hit.t_max == INFINITY)
            return;
        if (!TransparentTest(ray, hit)) {
            depth = hit.t_max;
            face = hit.face;
            return;
        }
        hit = HitRecord(hit.t_max + eps_zero, INFINITY);
    }
}