#include <iostream>
#include "model.h"

HitInfo::HitInfo() :
    position(vec3(NAN)), shapeNormal(NAN), surfaceNormal(NAN),
    opacity(1.0f), specular(0.04f), roughness(0.8f), metallic(0.0f), eta(1.0f),
    emission(0.0f), baseColor(0.0f), entering(true), id(0) {}

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
        if (Clum < eps_zero)
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

void Model::processMesh(aiMesh *mesh) {

    unsigned int offset = vertexDatas.size();
    for (unsigned int j = 0; j < mesh->mNumVertices; j++) {
        aiVector3D texCoord = mesh->mTextureCoords[0][j];
        aiVector3D normal = mesh->mNormals[j];
        vertexDatas.emplace_back(
                vec2(texCoord.x, texCoord.y),
                vec3(normal.x, normal.y, normal.z)
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
            transformedVertex[k] = vec3(vertex.x, vertex.y, vertex.z);
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
    }

    Face *meshFaces = &faces[offset];
    if (!material.texture[aiTextureType_EMISSIVE].empty())
        checkLightObject(meshFaces, mesh, material);
}

void Model::processMaterial(const std::string &model_folder, const aiScene *scene) {

    materials.resize(scene->mNumMaterials);
    for (int i = 0; i < scene->mNumMaterials; i++) {
        std::cout << "Loading material " << i << std::endl;
        aiMaterial *material = scene->mMaterials[i];

        aiString matName;
        if (material->Get(AI_MATKEY_NAME, matName) == AI_SUCCESS) {
            std::string name = matName.C_Str();
            std::cout << "Material Name: " << name << std::endl;
        }

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
        materials[i].id = i;
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
    importer.SetPropertyInteger(AI_CONFIG_PP_PTV_KEEP_HIERARCHY, 1);
    model_path = model_folder + model_name;
    const aiScene* scene = importer.ReadFile(model_path.c_str(),
                                             aiProcess_Triangulate |
                                             aiProcess_PreTransformVertices |
                                             aiProcess_SortByPType |
                                             aiProcess_FixInfacingNormals);

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

    for (int i = 0; i < scene->mNumMeshes; i++)
        processMesh(scene->mMeshes[i]);

    importer.FreeScene();

    kdt.build(faces);
}

bool TransparentTest(const Ray &ray, const HitRecord &hit) {
    const Face& face = *hit.face;
    const Material& material = *face.material;
    if (!material.hasFullyTransparentPart)
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

bool reverseFix(vec3 &v, const vec3 &Dir) {
    if (dot(v, Dir) < 0.0f) {
        v *= -1.0f;
        return false;
    }
    return true;
}

void getHitNormals(const Face& face, const vec3 &inDir, const vec3 &baryCoords,
                   vec3 &shapeNormal, vec3 &surfaceNormal_raw, bool &entering) {

    vec3 crossV0 = cross(face.v[1] - face.v[0], face.v[2] - face.v[0]);
    shapeNormal = normalize(crossV0);

    entering = reverseFix(shapeNormal, -inDir);
    surfaceNormal_raw = shapeNormal;

    float area = length(crossV0) / 2.0f;
    static const float areaThreshold = 1e-2f;
    if (area > areaThreshold)
        return;

    vec3 normalV0 = face.data[0]->normal;
    vec3 normalV1 = face.data[1]->normal;
    vec3 normalV2 = face.data[2]->normal;
    reverseFix(normalV0, shapeNormal);
    reverseFix(normalV1, shapeNormal);
    reverseFix(normalV2, shapeNormal);

    static const float threshold = 0.85f;
    surfaceNormal_raw = normalize(
            baryCoords[0] * (dot(normalV0, shapeNormal) > threshold ? normalV0 : shapeNormal)
            + baryCoords[1] * (dot(normalV1, shapeNormal) > threshold ? normalV1 : shapeNormal)
            + baryCoords[2] * (dot(normalV2, shapeNormal) > threshold ? normalV2 : shapeNormal)
    );

    if (!isfinite(surfaceNormal_raw)) {
        surfaceNormal_raw = shapeNormal;
        std::cerr << "  getHitNormals : NAN" << std::endl;
    }
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

void calcSurfaceNormal(const Face& face, const vec3 &normalMap, const vec3 &shapeNormal,
                       vec3 &surfaceNormal) {
    vec3 edge1 = face.v[1] - face.v[0];
    vec3 edge2 = face.v[2] - face.v[0];
    vec2 deltaUV1 = face.data[1]->uv - face.data[0]->uv;
    vec2 deltaUV2 = face.data[2]->uv - face.data[0]->uv;
    float f = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);
    vec3 tangentBasisU = f * (deltaUV2.y * edge1 - deltaUV1.y * edge2);
    vec3 tangentBasisV = f * (-deltaUV2.x * edge1 + deltaUV1.x * edge2);
    vec3 tangent = normalize(tangentBasisU - shapeNormal * dot(shapeNormal, tangentBasisU));
    vec3 bitangent = normalize(tangentBasisV - shapeNormal * dot(shapeNormal, tangentBasisV) - tangent * dot(tangent, tangentBasisV));
    vec3 sav = surfaceNormal;
    surfaceNormal = normalize(
            tangent * normalMap.x +
            bitangent * normalMap.y +
            surfaceNormal
    );
    if (!isfinite(surfaceNormal))
        surfaceNormal = sav;
}

void getHitMaterial(const Face& face, const vec3 &baryCoords, const vec3 &hit_dPdx, const vec3 &hit_dPdy,
                    HitInfo &hitInfo) {

    vec2 texUV =
            baryCoords[0] * face.data[0]->uv +
            baryCoords[1] * face.data[1]->uv +
            baryCoords[2] * face.data[2]->uv;

    const Material& material = *face.material;

    hitInfo.id = material.id;
    material.getSurfaceData(texUV[0], texUV[1], hitInfo.roughness, hitInfo.metallic);
    hitInfo.opacity = material.opacity;
    hitInfo.eta = material.ior; // 暂存绝对折射率，后续需要转换为相对折射率
    if (hitInfo.opacity > 1.0f - eps_zero)
        hitInfo.entering = true;
    vec2 dUVdx = getDuv(face, hit_dPdx),
         dUVdy = getDuv(face, hit_dPdy);
    float duv = (isfinite(dUVdx) && isfinite(dUVdx)) ?
            std::max((length(dUVdx)+length(dUVdy))/2.0f, sqrt(length(dUVdx)*length(dUVdy)))
            : NAN;

    if (hitInfo.opacity < eps_zero)
        hitInfo.baseColor = material.transmittingColor;
    else
        hitInfo.baseColor = material.getDiffuseColor(texUV[0], texUV[1], duv);
    hitInfo.emission = material.getEmissiveColor(texUV[0], texUV[1], duv);
    vec3 normalMap = material.getNormal(texUV[0], texUV[1], duv);
    calcSurfaceNormal(face, normalMap, hitInfo.shapeNormal, hitInfo.surfaceNormal);
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