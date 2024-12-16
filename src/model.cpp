#include <iostream>
#include "model.h"

HitInfo::HitInfo() {
    position = vec3(NAN);
}

unsigned int vertexCount(const aiScene *scene) {
    unsigned int cnt = 0;
    for (int i = 0; i < scene->mNumMeshes; i++)
        cnt += scene->mMeshes[i]->mNumVertices;
    return cnt;
}

unsigned int faceCount(const aiScene *scene) {
    unsigned int cnt = 0;
    for (int i = 0; i < scene->mNumMeshes; i++)
        cnt += scene->mMeshes[i]->mNumFaces;
    return cnt;
}

glm::mat4 aiMatrix4x4ToGlm(const aiMatrix4x4 &from) {
    glm::mat4 to;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            to[i][j] = from[j][i];
    return to;
}

vec3 getAverageEmissiveColor(const Material &material, const Face &face) {
    vec3 averageColor(0.0f);

    static const int gridResolution = 8;
    for (int i = 0; i < gridResolution; ++i) {
        for (int j = 0; j < gridResolution; ++j) {
            float a = static_cast<float>(i) / (gridResolution - 1);
            float b = static_cast<float>(j) / (gridResolution - 1);
            if (a + b > 1.0f) {
                a = 1.0f - a;
                b = 1.0f - b;
            }
            float c = 1.0f - a - b;
            vec2 uv = a * face.data[0]->uv + b * face.data[1]->uv + c * face.data[2]->uv;
            vec3 textureColor = material.getEmissiveColor(uv[0], uv[1], NAN);
            averageColor += textureColor;
        }
    }
    return averageColor / static_cast<float>(gridResolution * gridResolution);
}

void Model::checkLightObject(Face *meshFaces, aiMesh *mesh, const Material &material) {
    lightObjects.emplace_back();
    auto &lightObject = lightObjects.back();

    vec3 color = vec3(0.0f);

    std::vector<float> faceWeights;
    for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
        auto &face = meshFaces[j];

        vec3 textureColor = getAverageEmissiveColor(material, face);
        float Clum = dot(textureColor, RGB_Weight);
        if (Clum == 0.0f)
            continue;

        float area = length(cross(face.v[1] - face.v[0], face.v[2] - face.v[0])) / 2.0f;
        if (area < eps_zero)
            continue;

        color += textureColor * area;
        float power = area * Clum;
        lightObject.faces.emplace_back(face);
        faceWeights.emplace_back(power);
    }
    if (lightObject.faces.empty()) {
        lightObjects.pop_back();
        return ;
    }
    lightObject.powerDensity = 1.0f;
    lightObject.color = color / dot(color, RGB_Weight);
    for (unsigned int j = 0; j < lightObject.faces.size(); j++) {
        faceWeights[j] *= lightObject.powerDensity;
        lightObject.power += faceWeights[j];
        lightObject.center += lightObject.faces[j].center() * faceWeights[j];
    }

    lightObject.faceDist.Init(faceWeights);
    lightObject.center /= lightObject.power;

    std::cout << "Light object with power: " << lightObject.power
              << ", power density: " << lightObject.powerDensity
              << ", color:" << lightObject.color.x << ", " << lightObject.color.y << ", " << lightObject.color.z
              << ", position:" << lightObject.center.x << ", " << lightObject.center.y << ", "
              << lightObject.center.z << std::endl;
}

void Model::processMesh(aiMesh *mesh, const glm::mat4 &nodeTransform) {

    unsigned int offset = vertexDatas.size();
    for (unsigned int j = 0; j < mesh->mNumVertices; j++) {
        aiVector3D texCoord = mesh->mTextureCoords[0][j];
        vertexDatas.emplace_back(
                vec2(texCoord.x, texCoord.y),
                vec3(0)
        );
    }

    const auto &material = materials[mesh->mMaterialIndex];
    VertexData *vData = &vertexDatas[offset];
    offset = faces.size();
    for (unsigned int j = 0; j < mesh->mNumFaces; j++) {
        aiFace &face0 = mesh->mFaces[j];

        // Apply transformation to vertices
        glm::vec3 transformedVertex[3];
        for (int k = 0; k < 3; ++k) {
            aiVector3D vertex = mesh->mVertices[face0.mIndices[k]];
            transformedVertex[k] = (vec3)(nodeTransform * glm::vec4(vertex.x, vertex.y, vertex.z, 1.0f));
        }

        faces.push_back({
                {transformedVertex[0],
                 transformedVertex[1],
                 transformedVertex[2]},
                {&vData[face0.mIndices[0]],
                 &vData[face0.mIndices[1]],
                 &vData[face0.mIndices[2]]},
                &material,
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
    if (!material.texture[aiTextureType_EMISSIVE].empty())
        checkLightObject(meshFaces, mesh, material);
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

    processMaterial(model_folder, scene);

    unsigned int vertexCnt = vertexCount(scene);
    std::cout << "Vertices: " << vertexCnt << std::endl;
    vertexDatas.reserve(vertexCnt);
    unsigned int faceCnt = faceCount(scene);
    std::cout << "faces: " << faceCnt << std::endl;
    faces.reserve(faceCnt);

    const glm::mat4 identity = glm::mat4(1.0f);
    processNode(scene->mRootNode, scene, identity);

    importer.FreeScene();

    kdt.build(faces);
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
    vec4 diffuseColor = material.getDiffuseColor(texUV[0], texUV[1], NAN);
    return diffuseColor[3] < eps_zero;
}

void reverseFix(vec3 &v, const vec3 &Dir) {
    if (dot(v, Dir) < 0.0f)
        v *= -1.0f;
}

void getHitAllNormals(const Face& face, const vec3 &inDir, const vec3 &baryCoords,
                      vec3 &shapeNormal_to_compute, vec3 &surfaceNormal_to_compute_raw) {

    shapeNormal_to_compute = normalize(cross(face.v[1] - face.v[0], face.v[2] - face.v[0]));

    const auto &shapeNormal = shapeNormal_to_compute;

    reverseFix(shapeNormal_to_compute, -inDir);
    surfaceNormal_to_compute_raw = shapeNormal_to_compute;

    vec3 normalV0 = face.data[0]->normal;
    vec3 normalV1 = face.data[1]->normal;
    vec3 normalV2 = face.data[2]->normal;

    reverseFix(normalV0, shapeNormal_to_compute);
    reverseFix(normalV1, shapeNormal_to_compute);
    reverseFix(normalV2, shapeNormal_to_compute);
    static const float threshold = 0.9f;
    surfaceNormal_to_compute_raw = normalize(
            baryCoords[0] * (dot(normalV0, shapeNormal_to_compute) > threshold ? normalV0 : shapeNormal_to_compute)
            + baryCoords[1] * (dot(normalV1, shapeNormal_to_compute) > threshold ? normalV1 : shapeNormal_to_compute)
            + baryCoords[2] * (dot(normalV2, shapeNormal_to_compute) > threshold ? normalV2 : shapeNormal_to_compute)
    );
    if (isnan(surfaceNormal_to_compute_raw))
        surfaceNormal_to_compute_raw = shapeNormal_to_compute;
}

vec2 getDuv(const Face &face, const vec3 &hit_dPdx) {
    vec3 baryCoords = barycentric(face.v[0], face.v[1], face.v[2], face.v[0] + hit_dPdx);
    vec2 texUV0 = face.data[0]->uv;
    vec2 texUV1 =
            baryCoords[0] * face.data[0]->uv +
            baryCoords[1] * face.data[1]->uv +
            baryCoords[2] * face.data[2]->uv;
    return texUV1 - texUV0;
}

void getHitTexture(const Face& face, const vec3 &baryCoords, const vec3 &hit_dPdx, const vec3 &hit_dPdy, HitInfo &hitInfo_to_init) {

    // 插值得到命中点处的UV坐标
    vec2 texUV =
            baryCoords[0] * face.data[0]->uv +
            baryCoords[1] * face.data[1]->uv +
            baryCoords[2] * face.data[2]->uv;

    const Material& material = *face.material;

    vec3 &shapeNormal = hitInfo_to_init.shapeNormal,
         &surfaceNormal_to_compute = hitInfo_to_init.surfaceNormal;

    material.getSurfaceData(texUV[0], texUV[1], hitInfo_to_init.roughness, hitInfo_to_init.metallic);

    vec2 dUVdx = getDuv(face, hit_dPdx),
         dUVdy = getDuv(face, hit_dPdy);

    if (isnan(dUVdx) || isnan(dUVdx)) {
        surfaceNormal_to_compute = shapeNormal;
        vec3 normalMap = face.material->getNormal(texUV[0], texUV[1], NAN);
        hitInfo_to_init.baseColor = material.getDiffuseColor(texUV[0], texUV[1], NAN);
        hitInfo_to_init.emission = (material.texture[aiTextureType_EMISSIVE].empty())
                ? vec3(0.0f)
                : material.getEmissiveColor(texUV[0], texUV[1], NAN);
    } else {
        float duv = std::max(length(dUVdx), length(dUVdy));

        vec3 normalMap = face.material->getNormal(texUV[0], texUV[1], duv);
        hitInfo_to_init.baseColor = material.getDiffuseColor(texUV[0], texUV[1], duv);
        hitInfo_to_init.emission = (material.texture[aiTextureType_EMISSIVE].empty())
                ? vec3(0.0f)
                : material.getEmissiveColor(texUV[0], texUV[1], duv);

        // 根据UV差分求出边缘向量与UV差分
        vec3 edge1 = face.v[1] - face.v[0];
        vec3 edge2 = face.v[2] - face.v[0];
        vec2 deltaUV1 = face.data[1]->uv - face.data[0]->uv;
        vec2 deltaUV2 = face.data[2]->uv - face.data[0]->uv;
        float f = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);
        vec3 tangentBasisU = f * (deltaUV2.y * edge1 - deltaUV1.y * edge2);
        vec3 tangentBasisV = f * (-deltaUV2.x * edge1 + deltaUV1.x * edge2);
        vec3 tangent = normalize(tangentBasisU - shapeNormal * dot(shapeNormal, tangentBasisU));
        vec3 bitangent = normalize(tangentBasisV - shapeNormal * dot(shapeNormal, tangentBasisV) - tangent * dot(tangent, tangentBasisV));

        // 利用法线贴图，将TBN空间的法线转换到世界坐标系
        surfaceNormal_to_compute = normalize(
                tangent * normalMap.x +
                bitangent * normalMap.y +
                surfaceNormal_to_compute
        );

        if (isnan(surfaceNormal_to_compute))
            surfaceNormal_to_compute = shapeNormal;
    }
}

const int maxRayDepth_hit = 8;

HitRecord Model::rayHit(Ray ray) const {
    HitRecord hit(eps_zero, INFINITY);
    for (int T = 0; T < maxRayDepth_hit; T++) {
        kdt.rayHit(ray, hit);
        if (hit.t_max == INFINITY || !TransparentTest(ray, hit))
            return hit;
        hit = HitRecord(hit.t_max + eps_zero, INFINITY);
    }
    return hit;
}

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