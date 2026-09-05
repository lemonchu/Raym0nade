#include "raym0nade/scene_data.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "raym0nade/model.hpp"

namespace raym0nade {
namespace {

template <typename Value>
std::size_t borrowedIndex(
    const Value* pointer, const std::vector<Value>& values, const char* description) {
    if (pointer == nullptr || values.empty()) {
        throw std::logic_error(std::string{"Packed scene contains a missing "} + description + '.');
    }
    const Value* begin = values.data();
    const Value* end = begin + values.size();
    const std::less<const Value*> less;
    if (less(pointer, begin) || !less(pointer, end)) {
        throw std::logic_error(
            std::string{"Packed scene contains an out-of-range "} + description + '.');
    }
    return static_cast<std::size_t>(pointer - begin);
}

bool allFinite(const std::array<float, 4>& values) noexcept {
    return std::all_of(values.begin(), values.end(), [](float value) {
        return std::isfinite(value);
    });
}

bool samePosition(const PackedVertex& vertex, const vec3& position) noexcept {
    return vertex.positionAndNormalX[0] == position.x &&
           vertex.positionAndNormalX[1] == position.y &&
           vertex.positionAndNormalX[2] == position.z;
}

PackedVertex packVertex(const vec3& position, const VertexData& attributes) {
    if (!isFinite(position) || !isFinite(attributes.normal) || !isFinite(attributes.uv)) {
        throw std::invalid_argument("Packed scene vertices must contain only finite values.");
    }
    return PackedVertex{
        {position.x, position.y, position.z, attributes.normal.x},
        {attributes.normal.y, attributes.normal.z, attributes.uv.x, attributes.uv.y},
    };
}

PackedMaterial packMaterial(const Material& material) {
    PackedMaterial result;
    result.diffuseAndOpacity = {
        material.diffuseFactor.x,
        material.diffuseFactor.y,
        material.diffuseFactor.z,
        material.opacity,
    };
    result.emissionAndIor = {
        material.emissiveFactor.x,
        material.emissiveFactor.y,
        material.emissiveFactor.z,
        material.ior,
    };
    result.transmissionAndRoughness = {
        material.transmissionColor.x,
        material.transmissionColor.y,
        material.transmissionColor.z,
        material.roughness,
    };
    result.metallicSpecularAndReserved = {material.metallic, 0.04F, 0.0F, 0.0F};
    if (material.hasCutoutTransparency) {
        result.flagsAndReserved[0] |= kPackedMaterialCutout;
    }
    if (material.hasTexture(TextureSlot::diffuse)) {
        result.flagsAndReserved[0] |= kPackedMaterialHasDiffuseTexture;
    }
    if (material.hasTexture(TextureSlot::specular)) {
        result.flagsAndReserved[0] |= kPackedMaterialHasSpecularTexture;
    }
    if (material.hasTexture(TextureSlot::emissive)) {
        result.flagsAndReserved[0] |= kPackedMaterialHasEmissiveTexture;
    }
    if (material.hasTexture(TextureSlot::normal)) {
        result.flagsAndReserved[0] |= kPackedMaterialHasNormalTexture;
    }
    return result;
}

}  // namespace

std::size_t PackedSceneData::triangleCount() const noexcept {
    return triangleIndices.size() / 3U;
}

void PackedSceneData::validate() const {
    if (formatVersion != kPackedSceneFormatVersion) {
        throw std::invalid_argument("Packed scene format version is unsupported.");
    }
    if (vertices.empty() || triangleIndices.empty() || materials.empty()) {
        throw std::invalid_argument("Packed scene geometry and materials must not be empty.");
    }
    if (triangleIndices.size() % 3U != 0U ||
        triangleMaterialIds.size() != triangleCount()) {
        throw std::invalid_argument("Packed scene triangle arrays have inconsistent sizes.");
    }
    if (vertices.size() > kInvalidSceneId || materials.size() > kInvalidSceneId ||
        triangleCount() > kInvalidSceneId) {
        throw std::invalid_argument("Packed scene counts exceed the 32-bit index range.");
    }
    for (const PackedVertex& vertex : vertices) {
        if (!allFinite(vertex.positionAndNormalX) || !allFinite(vertex.normalYZAndUv)) {
            throw std::invalid_argument("Packed scene vertices must contain only finite values.");
        }
    }
    for (std::uint32_t index : triangleIndices) {
        if (index >= vertices.size()) {
            throw std::invalid_argument("Packed scene contains an out-of-range vertex index.");
        }
    }
    for (std::uint32_t materialId : triangleMaterialIds) {
        if (materialId >= materials.size()) {
            throw std::invalid_argument("Packed scene contains an out-of-range material index.");
        }
    }
    for (const PackedMaterial& material : materials) {
        if (!allFinite(material.diffuseAndOpacity) || !allFinite(material.emissionAndIor) ||
            !allFinite(material.transmissionAndRoughness) ||
            !allFinite(material.metallicSpecularAndReserved)) {
            throw std::invalid_argument("Packed scene materials must contain only finite values.");
        }
        if ((material.flagsAndReserved[0] & ~kPackedMaterialKnownFlags) != 0U) {
            throw std::invalid_argument("Packed scene material contains unknown flags.");
        }
        if (material.flagsAndReserved[1] != 0U || material.flagsAndReserved[2] != 0U ||
            material.flagsAndReserved[3] != 0U ||
            material.metallicSpecularAndReserved[2] != 0.0F ||
            material.metallicSpecularAndReserved[3] != 0.0F) {
            throw std::invalid_argument("Packed scene material reserved fields must be zero.");
        }
        if (std::any_of(
                material.textureIds.begin(),
                material.textureIds.end(),
                [](std::uint32_t textureId) { return textureId != kInvalidSceneId; })) {
            throw std::invalid_argument(
                "Packed scene texture IDs must remain invalid until texture packing is supported.");
        }
    }
}

PackedSceneData Model::packScene() const {
    if (vertexData_.size() >= kInvalidSceneId || materials_.size() >= kInvalidSceneId ||
        faces_.size() >= kInvalidSceneId ||
        faces_.size() > std::numeric_limits<std::size_t>::max() / 3U) {
        throw std::overflow_error("Model exceeds the packed scene index range.");
    }

    PackedSceneData result;
    result.vertices.reserve(vertexData_.size());
    result.triangleIndices.reserve(faces_.size() * 3U);
    result.triangleMaterialIds.reserve(faces_.size());
    result.materials.reserve(materials_.size());
    for (const Material& material : materials_) {
        result.materials.push_back(packMaterial(material));
    }

    std::vector<std::uint32_t> sourceToPacked(vertexData_.size(), kInvalidSceneId);
    for (const Face& face : faces_) {
        const std::size_t materialIndex = borrowedIndex(face.material, materials_, "material pointer");
        result.triangleMaterialIds.push_back(static_cast<std::uint32_t>(materialIndex));
        for (std::size_t corner = 0; corner < 3U; ++corner) {
            const std::size_t sourceIndex =
                borrowedIndex(face.vertexData[corner], vertexData_, "vertex pointer");
            std::uint32_t& packedIndex = sourceToPacked[sourceIndex];
            if (packedIndex == kInvalidSceneId) {
                if (result.vertices.size() >= kInvalidSceneId) {
                    throw std::overflow_error("Packed scene has too many referenced vertices.");
                }
                packedIndex = static_cast<std::uint32_t>(result.vertices.size());
                result.vertices.push_back(packVertex(face.vertices[corner], vertexData_[sourceIndex]));
            } else if (!samePosition(result.vertices[packedIndex], face.vertices[corner])) {
                throw std::logic_error(
                    "A source vertex maps to inconsistent positions in the packed scene.");
            }
            result.triangleIndices.push_back(packedIndex);
        }
    }

    result.validate();
    return result;
}

}  // namespace raym0nade
