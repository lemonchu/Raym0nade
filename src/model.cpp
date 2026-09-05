#include "raym0nade/model.hpp"

#include <assimp/Importer.hpp>
#include <assimp/config.h>
#include <assimp/material.h>
#include <assimp/postprocess.h>
#include <assimp/scene.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace raym0nade {
namespace {

constexpr int kEmissiveSampleGrid = 8;
constexpr int kMaximumCutoutLayers = 32;
constexpr float kMinimumEmitterImportance = 1.0e-6F;

std::uint64_t nextModelInstanceIdentity() {
    static std::atomic<std::uint64_t> nextIdentity{1U};
    std::uint64_t current = nextIdentity.load(std::memory_order_relaxed);
    while (current != std::numeric_limits<std::uint64_t>::max()) {
        if (nextIdentity.compare_exchange_weak(
                current,
                current + 1U,
                std::memory_order_relaxed,
                std::memory_order_relaxed)) {
            return current;
        }
    }
    throw std::overflow_error(std::string{});
}

std::size_t countSceneVertices(const aiScene& scene) noexcept {
    std::size_t count = 0;
    for (unsigned int index = 0; index < scene.mNumMeshes; ++index) {
        count += scene.mMeshes[index]->mNumVertices;
    }
    return count;
}

std::size_t countSceneFaces(const aiScene& scene) noexcept {
    std::size_t count = 0;
    for (unsigned int index = 0; index < scene.mNumMeshes; ++index) {
        count += scene.mMeshes[index]->mNumFaces;
    }
    return count;
}

vec3 toVector(const aiVector3D& value) noexcept {
    return vec3{value.x, value.y, value.z};
}

vec3 toColor(const aiColor3D& value, const vec3& fallback) noexcept {
    const vec3 color{value.r, value.g, value.b};
    return isFinite(color) ? glm::max(color, vec3{0.0F}) : fallback;
}

std::filesystem::path texturePath(
    const std::filesystem::path& modelDirectory, const std::string& encodedTexturePath) {
    std::string decoded = urlDecode(encodedTexturePath);
    std::replace(decoded.begin(), decoded.end(), '\\', '/');
    const std::filesystem::path path{decoded};
    return (path.is_absolute() ? path : modelDirectory / path).lexically_normal();
}

vec3 averageEmissiveColor(const Material& material, const Face& face) {
    vec3 average{0.0F};
    for (int row = 0; row < kEmissiveSampleGrid; ++row) {
        for (int column = 0; column < kEmissiveSampleGrid; ++column) {
            float a = static_cast<float>(row) / static_cast<float>(kEmissiveSampleGrid - 1);
            float b = static_cast<float>(column) / static_cast<float>(kEmissiveSampleGrid - 1);
            if (a + b > 1.0F) {
                a = 1.0F - a;
                b = 1.0F - b;
            }
            const float c = 1.0F - a - b;
            const vec2 uv = a * face.vertexData[0]->uv + b * face.vertexData[1]->uv +
                            c * face.vertexData[2]->uv;
            average += material.emissiveColor(uv.x, uv.y, std::numeric_limits<float>::quiet_NaN());
        }
    }
    return average / static_cast<float>(kEmissiveSampleGrid * kEmissiveSampleGrid);
}

bool orientToward(vec3& normal, const vec3& direction) noexcept {
    if (glm::dot(normal, direction) < 0.0F) {
        normal = -normal;
        return false;
    }
    return true;
}

bool isTransparentCutout(const Ray& ray, const HitRecord& hit) {
    if (hit.face == nullptr || hit.face->material == nullptr ||
        !hit.face->material->hasCutoutTransparency) {
        return false;
    }

    const Face& face = *hit.face;
    const vec3 position = ray.origin + ray.direction * hit.tMaximum;
    const vec3 coordinates =
        barycentric(face.vertices[0], face.vertices[1], face.vertices[2], position);
    if (!isFinite(coordinates)) {
        return false;
    }
    const vec2 uv = coordinates[0] * face.vertexData[0]->uv +
                    coordinates[1] * face.vertexData[1]->uv +
                    coordinates[2] * face.vertexData[2]->uv;
    return face.material->diffuseColor(
               uv.x, uv.y, std::numeric_limits<float>::quiet_NaN()).a < kRayEpsilon;
}

vec2 textureDerivative(const Face& face, const vec3& positionDelta) noexcept {
    const vec3 coordinates = barycentric(
        face.vertices[0], face.vertices[1], face.vertices[2], face.vertices[0] + positionDelta);
    if (!isFinite(coordinates)) {
        return vec2{std::numeric_limits<float>::quiet_NaN()};
    }
    return (coordinates[0] - 1.0F) * face.vertexData[0]->uv +
           coordinates[1] * face.vertexData[1]->uv + coordinates[2] * face.vertexData[2]->uv;
}

void applyNormalMap(
    const Face& face, const vec3& normalMap, const vec3& shapeNormal, vec3& surfaceNormal) noexcept {
    if (!isFinite(normalMap) || glm::dot(normalMap, normalMap) <= kRayEpsilon * kRayEpsilon) {
        return;
    }

    const vec3 edge1 = face.vertices[1] - face.vertices[0];
    const vec3 edge2 = face.vertices[2] - face.vertices[0];
    const vec2 deltaUv1 = face.vertexData[1]->uv - face.vertexData[0]->uv;
    const vec2 deltaUv2 = face.vertexData[2]->uv - face.vertexData[0]->uv;
    const float determinant = deltaUv1.x * deltaUv2.y - deltaUv2.x * deltaUv1.y;
    if (!isFinite(determinant) || std::abs(determinant) <= 1.0e-8F) {
        return;
    }

    const float inverseDeterminant = 1.0F / determinant;
    const vec3 tangentCandidate = inverseDeterminant * (deltaUv2.y * edge1 - deltaUv1.y * edge2);
    const vec3 tangent = safeNormalize(
        tangentCandidate - shapeNormal * glm::dot(shapeNormal, tangentCandidate));
    if (glm::dot(tangent, tangent) <= 0.0F) {
        return;
    }
    const vec3 bitangentCandidate = inverseDeterminant *
                                    (-deltaUv2.x * edge1 + deltaUv1.x * edge2);
    const float handedness = glm::dot(glm::cross(shapeNormal, tangent), bitangentCandidate) < 0.0F
                                 ? -1.0F
                                 : 1.0F;
    const vec3 bitangent = handedness * safeNormalize(glm::cross(shapeNormal, tangent));
    const vec3 mapped = tangent * normalMap.x + bitangent * normalMap.y + surfaceNormal * normalMap.z;
    const vec3 normalized = safeNormalize(mapped, surfaceNormal);
    if (isFinite(normalized)) {
        surfaceNormal = normalized;
    }
}

}  // namespace

class ModelBuilder {
public:
    explicit ModelBuilder(Model& model) noexcept : model_(model) {}

    void build();

private:
    static void loadMaterialProperties(Material& destination, const aiMaterial& source);
    void processMaterials(const aiScene& scene);
    void processMesh(const aiMesh& mesh);
    void createLightObject(Face* meshFaces, std::size_t meshFaceCount, const Material& material);

    Model& model_;
};

Model::Model(
    const std::filesystem::path& modelDirectory,
    const std::filesystem::path& modelFilename,
    const std::filesystem::path& skyFilename)
    : instanceIdentity_(nextModelInstanceIdentity()) {
    if (modelDirectory.empty() || modelFilename.empty()) {
        throw std::invalid_argument("Model directory and filename must not be empty.");
    }

    modelPath_ =
        (modelFilename.is_absolute() ? modelFilename : modelDirectory / modelFilename)
            .lexically_normal();
    if (!skyFilename.empty() && skyFilename != "null") {
        skyPath_ = (skyFilename.is_absolute() ? skyFilename
                                              : modelPath_.parent_path() / skyFilename)
                       .lexically_normal();
        std::cout << "Loading sky map: " << skyPath_.u8string() << '\n';
        sky_.load(skyPath_);
    }

    ModelBuilder{*this}.build();
}

void ModelBuilder::build() {
    Assimp::Importer importer;
    importer.SetPropertyInteger(AI_CONFIG_PP_PTV_KEEP_HIERARCHY, 1);
    // The first FBX geometry layer contains the Bistro renderable geometry and avoids
    // importing duplicate layers. Assimp 5.4 can still omit meshes nondeterministically;
    // appearance and benchmark runs record their imported topology instead of hiding it.
    importer.SetPropertyBool(AI_CONFIG_IMPORT_FBX_READ_ALL_GEOMETRY_LAYERS, false);
    const unsigned int flags = aiProcess_Triangulate | aiProcess_PreTransformVertices;
    const std::string modelPathUtf8 = model_.modelPath_.u8string();
    const aiScene* scene = importer.ReadFile(modelPathUtf8, flags);
    if (scene == nullptr || (scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE) != 0U || scene->mRootNode == nullptr) {
        throw std::runtime_error(
            "Failed to load model " + modelPathUtf8 + ": " + importer.GetErrorString());
    }

    const unsigned int sceneMeshCount = scene->mNumMeshes;
    const std::size_t sceneVertexCount = countSceneVertices(*scene);
    const std::size_t sceneFaceCount = countSceneFaces(*scene);
    processMaterials(*scene);
    if (scene->mNumMeshes != sceneMeshCount || countSceneVertices(*scene) != sceneVertexCount ||
        countSceneFaces(*scene) != sceneFaceCount) {
        throw std::runtime_error("The imported scene topology changed while loading materials.");
    }
    model_.vertexData_.reserve(sceneVertexCount);
    model_.faces_.reserve(sceneFaceCount);
    std::cout << "Vertices: " << sceneVertexCount << '\n';
    std::cout << "Faces: " << sceneFaceCount << '\n';
    for (unsigned int index = 0; index < scene->mNumMeshes; ++index) {
        processMesh(*scene->mMeshes[index]);
    }

    if (model_.faces_.empty()) {
        throw std::runtime_error("The model contains no renderable triangles: " + modelPathUtf8);
    }
    model_.bvh_.build(model_.faces_);
}

void ModelBuilder::loadMaterialProperties(Material& destination, const aiMaterial& source) {
    aiString materialName;
    if (source.Get(AI_MATKEY_NAME, materialName) == AI_SUCCESS) {
        destination.name = materialName.C_Str();
    }

    aiColor3D sourceDiffuse{1.0F, 1.0F, 1.0F};
    if (source.Get(AI_MATKEY_COLOR_DIFFUSE, sourceDiffuse) == AI_SUCCESS) {
        destination.diffuseFactor = toColor(sourceDiffuse, vec3{1.0F});
    }

    aiColor3D sourceEmission{0.0F, 0.0F, 0.0F};
    const vec3 importedEmission = source.Get(AI_MATKEY_COLOR_EMISSIVE, sourceEmission) == AI_SUCCESS
                                      ? toColor(sourceEmission, vec3{0.0F})
                                      : vec3{0.0F};
    destination.emissiveFactor =
        glm::dot(importedEmission, kLuminanceWeights) > 0.0F
            ? importedEmission
            : (destination.hasTexture(TextureSlot::emissive) ? vec3{1.0F} : vec3{0.0F});

    ai_real importedIor = destination.ior;
    if (source.Get(AI_MATKEY_REFRACTI, importedIor) == AI_SUCCESS &&
        isFinite(static_cast<float>(importedIor)) && importedIor > 0.0F) {
        destination.ior = static_cast<float>(importedIor);
    }

    ai_real importedRoughness = destination.roughness;
    if (source.Get(AI_MATKEY_ROUGHNESS_FACTOR, importedRoughness) == AI_SUCCESS &&
        isFinite(static_cast<float>(importedRoughness))) {
        destination.roughness =
            std::clamp(static_cast<float>(importedRoughness), 1.0e-3F, 1.0F);
    }

    ai_real importedMetallic = destination.metallic;
    if (source.Get(AI_MATKEY_METALLIC_FACTOR, importedMetallic) == AI_SUCCESS &&
        isFinite(static_cast<float>(importedMetallic))) {
        destination.metallic =
            std::clamp(static_cast<float>(importedMetallic), 0.0F, 1.0F);
    }

    ai_real sourceOpacity = 1.0F;
    if (source.Get(AI_MATKEY_OPACITY, sourceOpacity) != AI_SUCCESS ||
        !isFinite(static_cast<float>(sourceOpacity)) || sourceOpacity >= 0.99F) {
        return;
    }

    destination.opacity = 0.0F;
    if (destination.ior <= 1.0F + kRayEpsilon) {
        destination.ior = 1.25F;
    }
    destination.roughness = destination.name == "Ice" ? 0.5F : 5.0e-3F;
    destination.transmissionColor = destination.diffuseFactor;
    if (glm::dot(destination.transmissionColor, kLuminanceWeights) <= kRayEpsilon) {
        destination.transmissionColor = vec3{1.0F};
    }

    if (destination.name == "TransparentGlassWine") {
        destination.transmissionColor = vec3{0.2F, 0.08F, 0.07F};
    } else if (destination.name == "TransparentGlass" || destination.name == "Water" ||
               destination.name == "Ice") {
        destination.transmissionColor = vec3{1.0F};
    } else if (destination.name == "Beer") {
        destination.transmissionColor = vec3{0.8F, 0.7F, 0.55F};
    } else if (destination.name == "Red_Wine") {
        destination.transmissionColor = vec3{0.24F, 0.09F, 0.07F};
    } else if (destination.name == "White_Wine") {
        destination.transmissionColor = vec3{0.85F, 0.78F, 0.6F};
    }
}

void ModelBuilder::processMaterials(const aiScene& scene) {
    model_.materials_.resize(scene.mNumMaterials);
    const std::pair<aiTextureType, TextureSlot> textureTypes[] = {
        {aiTextureType_DIFFUSE, TextureSlot::diffuse},
        {aiTextureType_SPECULAR, TextureSlot::specular},
        {aiTextureType_EMISSIVE, TextureSlot::emissive},
        {aiTextureType_NORMALS, TextureSlot::normal},
    };

    for (unsigned int materialIndex = 0; materialIndex < scene.mNumMaterials; ++materialIndex) {
        Material& destination = model_.materials_[materialIndex];
        const aiMaterial& source = *scene.mMaterials[materialIndex];
        destination.id = static_cast<int>(materialIndex);

        for (const auto& [sourceType, destinationSlot] : textureTypes) {
            if (source.GetTextureCount(sourceType) == 0U) {
                continue;
            }
            aiString sourcePath;
            if (source.GetTexture(sourceType, 0, &sourcePath) != AI_SUCCESS) {
                continue;
            }
            const std::filesystem::path path =
                texturePath(model_.modelPath_.parent_path(), sourcePath.C_Str());
            std::cout << "Loading texture: " << path.u8string() << '\n';
            try {
                destination.loadTexture(destinationSlot, path);
            } catch (const std::exception& error) {
                std::cerr << "Texture skipped: " << error.what() << '\n';
            }
        }
        loadMaterialProperties(destination, source);
    }
    std::cout << "Materials: " << model_.materials_.size() << '\n';
}

void ModelBuilder::processMesh(const aiMesh& mesh) {
    if (mesh.mMaterialIndex >= model_.materials_.size()) {
        throw std::runtime_error("Mesh references an invalid material index.");
    }

    const std::size_t vertexOffset = model_.vertexData_.size();
    for (unsigned int index = 0; index < mesh.mNumVertices; ++index) {
        const vec2 uv = mesh.HasTextureCoords(0) ? vec2{mesh.mTextureCoords[0][index].x,
                                                       mesh.mTextureCoords[0][index].y}
                                                 : vec2{0.0F};
        const vec3 normal = mesh.HasNormals() ? safeNormalize(toVector(mesh.mNormals[index])) : vec3{0.0F};
        model_.vertexData_.push_back(VertexData{uv, normal});
    }

    Material& material = model_.materials_[mesh.mMaterialIndex];
    const std::size_t firstMeshFace = model_.faces_.size();
    std::size_t nonFiniteFaceCount = 0;
    for (unsigned int index = 0; index < mesh.mNumFaces; ++index) {
        const aiFace& sourceFace = mesh.mFaces[index];
        if (sourceFace.mNumIndices != 3U) {
            continue;
        }
        const unsigned int i0 = sourceFace.mIndices[0];
        const unsigned int i1 = sourceFace.mIndices[1];
        const unsigned int i2 = sourceFace.mIndices[2];
        if (i0 >= mesh.mNumVertices || i1 >= mesh.mNumVertices || i2 >= mesh.mNumVertices) {
            throw std::runtime_error("Mesh contains an out-of-range vertex index.");
        }
        const vec3 vertices[3]{
            toVector(mesh.mVertices[i0]),
            toVector(mesh.mVertices[i1]),
            toVector(mesh.mVertices[i2]),
        };
        if (!isFinite(vertices[0]) || !isFinite(vertices[1]) || !isFinite(vertices[2])) {
            ++nonFiniteFaceCount;
            continue;
        }
        model_.faces_.push_back(Face{
            {vertices[0], vertices[1], vertices[2]},
            {&model_.vertexData_[vertexOffset + i0],
             &model_.vertexData_[vertexOffset + i1],
             &model_.vertexData_[vertexOffset + i2]},
            &material,
        });
    }
    if (nonFiniteFaceCount > 0) {
        std::cerr << "Skipped " << nonFiniteFaceCount
                  << " triangle(s) with non-finite vertex coordinates.\n";
    }

    const std::size_t meshFaceCount = model_.faces_.size() - firstMeshFace;
    if (material.isEmissive() && meshFaceCount > 0) {
        createLightObject(&model_.faces_[firstMeshFace], meshFaceCount, material);
    }
}

void ModelBuilder::createLightObject(
    Face* meshFaces, std::size_t meshFaceCount, const Material& material) {
    LightObject light;
    vec3 weightedColor{0.0F};
    std::vector<float> faceWeights;
    faceWeights.reserve(meshFaceCount);

    for (std::size_t index = 0; index < meshFaceCount; ++index) {
        const Face& face = meshFaces[index];
        const vec3 textureColor = averageEmissiveColor(material, face);
        const float luminance = glm::dot(textureColor, kLuminanceWeights);
        const float area = 0.5F * glm::length(glm::cross(
                                      face.vertices[1] - face.vertices[0],
                                      face.vertices[2] - face.vertices[0]));
        if (!isFinite(luminance) || !isFinite(area) || area <= kRayEpsilon) {
            continue;
        }
        const float importance = material.hasTexture(TextureSlot::emissive)
                                     ? std::max(luminance, kMinimumEmitterImportance)
                                     : luminance;
        if (importance <= 0.0F) {
            continue;
        }
        const float power = area * importance;
        weightedColor += textureColor * area;
        light.faces.push_back(face);
        faceWeights.push_back(power);
    }

    if (light.faces.empty()) {
        return;
    }
    const float colorLuminance = glm::dot(weightedColor, kLuminanceWeights);
    light.color = colorLuminance > 0.0F && isFinite(colorLuminance)
                      ? weightedColor / colorLuminance
                      : vec3{1.0F};
    for (std::size_t index = 0; index < light.faces.size(); ++index) {
        light.power += faceWeights[index];
        light.center += light.faces[index].center() * faceWeights[index];
    }
    if (light.power <= 0.0F || !isFinite(light.power)) {
        return;
    }
    light.center /= light.power;
    light.faceDistribution.initialize(faceWeights);
    model_.lights_.push_back(std::move(light));
}

void getHitNormals(
    const Face& face,
    const vec3& incomingDirection,
    const vec3& barycentricCoordinates,
    vec3& shapeNormal,
    vec3& surfaceNormal,
    bool& entering) noexcept {
    const vec3 geometricNormal = glm::cross(
        face.vertices[1] - face.vertices[0], face.vertices[2] - face.vertices[0]);
    shapeNormal = safeNormalize(geometricNormal, -incomingDirection);
    entering = orientToward(shapeNormal, -incomingDirection);

    vec3 normals[] = {
        face.vertexData[0] != nullptr
            ? safeNormalize(face.vertexData[0]->normal, shapeNormal)
            : shapeNormal,
        face.vertexData[1] != nullptr
            ? safeNormalize(face.vertexData[1]->normal, shapeNormal)
            : shapeNormal,
        face.vertexData[2] != nullptr
            ? safeNormalize(face.vertexData[2]->normal, shapeNormal)
            : shapeNormal,
    };
    for (vec3& normal : normals) {
        orientToward(normal, shapeNormal);
    }
    const vec3 interpolated = barycentricCoordinates[0] * normals[0] +
                              barycentricCoordinates[1] * normals[1] +
                              barycentricCoordinates[2] * normals[2];
    surfaceNormal = safeNormalize(interpolated, shapeNormal);
}

void getHitMaterial(
    const Face& face,
    const vec3& barycentricCoordinates,
    const vec3& positionDx,
    const vec3& positionDy,
    HitInfo& hitInfo) {
    if (face.material == nullptr) {
        throw std::runtime_error("Hit face has no material.");
    }

    const vec2 uv = barycentricCoordinates[0] * face.vertexData[0]->uv +
                    barycentricCoordinates[1] * face.vertexData[1]->uv +
                    barycentricCoordinates[2] * face.vertexData[2]->uv;
    const Material& material = *face.material;
    hitInfo.materialId = material.id;
    material.surfaceParameters(uv.x, uv.y, hitInfo.roughness, hitInfo.metallic);
    hitInfo.opacity = material.opacity;
    hitInfo.eta = material.ior;
    if (hitInfo.opacity > 1.0F - kRayEpsilon) {
        hitInfo.entering = true;
    }

    const vec2 uvDx = textureDerivative(face, positionDx);
    const vec2 uvDy = textureDerivative(face, positionDy);
    const float footprint = isFinite(uvDx) && isFinite(uvDy)
                                ? 0.5F * (glm::length(uvDx) + glm::length(uvDy))
                                : std::numeric_limits<float>::quiet_NaN();
    hitInfo.baseColor = hitInfo.opacity < kRayEpsilon
                            ? material.transmissionColor
                            : vec3{material.diffuseColor(uv.x, uv.y, footprint)};
    hitInfo.emission = material.emissiveColor(uv.x, uv.y, footprint);
    applyNormalMap(face, material.normal(uv.x, uv.y, footprint), hitInfo.shapeNormal, hitInfo.surfaceNormal);
}

HitRecord Model::intersect(const Ray& ray) const noexcept {
    HitRecord hit{kRayEpsilon, std::numeric_limits<float>::infinity()};
    for (int layer = 0; layer < kMaximumCutoutLayers; ++layer) {
        bvh_.intersect(ray, hit);
        if (hit.face == nullptr || !isTransparentCutout(ray, hit)) {
            return hit;
        }
        hit = HitRecord{hit.tMaximum + kRayEpsilon, std::numeric_limits<float>::infinity()};
    }
    return HitRecord{kRayEpsilon, std::numeric_limits<float>::infinity()};
}

bool Model::occluded(const Ray& ray, float maximumDistance) const noexcept {
    if (maximumDistance <= kRayEpsilon) {
        return false;
    }
    HitRecord hit{kRayEpsilon, maximumDistance + kRayEpsilon};
    for (int layer = 0; layer < kMaximumCutoutLayers; ++layer) {
        bvh_.intersect(ray, hit);
        if (hit.face == nullptr || hit.tMaximum >= maximumDistance) {
            return false;
        }
        if (!isTransparentCutout(ray, hit)) {
            return true;
        }
        hit = HitRecord{hit.tMaximum + kRayEpsilon, maximumDistance + kRayEpsilon};
    }
    return false;
}

const std::vector<LightObject>& Model::lights() const noexcept {
    return lights_;
}

const SkyBox& Model::sky() const noexcept {
    return sky_;
}

const std::filesystem::path& Model::modelPath() const noexcept {
    return modelPath_;
}

std::uint64_t Model::instanceIdentity() const noexcept {
    return instanceIdentity_;
}

std::size_t Model::faceCount() const noexcept {
    return faces_.size();
}

}  // namespace raym0nade
