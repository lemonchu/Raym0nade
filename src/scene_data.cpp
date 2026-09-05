#include "raym0nade/scene_data.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
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

constexpr std::array<TextureSlot, 4> kTextureSlots{
    TextureSlot::diffuse,
    TextureSlot::specular,
    TextureSlot::emissive,
    TextureSlot::normal,
};

constexpr std::array<std::uint32_t, 4> kTextureFlags{
    kPackedMaterialHasDiffuseTexture,
    kPackedMaterialHasSpecularTexture,
    kPackedMaterialHasEmissiveTexture,
    kPackedMaterialHasNormalTexture,
};

static_assert(static_cast<std::size_t>(TextureSlot::diffuse) == 0U);
static_assert(static_cast<std::size_t>(TextureSlot::specular) == 1U);
static_assert(static_cast<std::size_t>(TextureSlot::emissive) == 2U);
static_assert(static_cast<std::size_t>(TextureSlot::normal) == 3U);
static_assert(static_cast<std::size_t>(TextureSlot::count) == kTextureSlots.size());

std::string normalizedTextureKey(const std::filesystem::path& sourcePath) {
    std::string key = sourcePath.lexically_normal().generic_u8string();
    if (key.empty()) {
        throw std::logic_error("A packed texture is missing its source path.");
    }
    return key;
}

std::uint8_t expandedComponent(
    const std::vector<std::uint8_t>& pixels,
    std::size_t pixelIndex,
    int channels,
    int requestedChannel) {
    const std::size_t base = pixelIndex * static_cast<std::size_t>(channels);
    if (channels == 1) {
        return requestedChannel == 3 ? 255U : pixels[base];
    }
    if (channels == 2) {
        return requestedChannel == 3 ? pixels[base + 1U] : pixels[base];
    }
    if (requestedChannel == 3) {
        return channels == 4 ? pixels[base + 3U] : 255U;
    }
    return pixels[base + static_cast<std::size_t>(requestedChannel)];
}

std::uint32_t packRgba8(
    const std::vector<std::uint8_t>& pixels, std::size_t pixelIndex, int channels) {
    const std::uint32_t red = expandedComponent(pixels, pixelIndex, channels, 0);
    const std::uint32_t green = expandedComponent(pixels, pixelIndex, channels, 1);
    const std::uint32_t blue = expandedComponent(pixels, pixelIndex, channels, 2);
    const std::uint32_t alpha = expandedComponent(pixels, pixelIndex, channels, 3);
    return red | (green << 8U) | (blue << 16U) | (alpha << 24U);
}

std::uint32_t appendTexture(PackedSceneData& scene, const ImageData& image) {
    if (image.empty() || image.width() <= 0 || image.height() <= 0 ||
        image.channels() < 1 || image.channels() > 4 || image.mipLevelCount() == 0U ||
        image.mipLevelCount() > kPackedTextureMaxMipLevels) {
        throw std::logic_error("A material contains invalid texture data.");
    }
    if (scene.textures.size() >= kInvalidSceneId) {
        throw std::overflow_error("Packed scene has too many textures.");
    }

    const std::size_t levelCount = image.mipLevelCount();
    const std::size_t maximumUint32 = std::numeric_limits<std::uint32_t>::max();
    if (scene.textureMipLevels.size() > maximumUint32 - levelCount) {
        throw std::overflow_error("Packed scene texture mip count exceeds the 32-bit range.");
    }

    std::size_t addedTexels = 0U;
    for (std::size_t level = 0U; level < levelCount; ++level) {
        const int mipWidth = image.mipWidth(level);
        const int mipHeight = image.mipHeight(level);
        if (mipWidth <= 0 || mipHeight <= 0) {
            throw std::logic_error("A material texture has an invalid mip extent.");
        }
        const std::size_t width = static_cast<std::size_t>(mipWidth);
        const std::size_t height = static_cast<std::size_t>(mipHeight);
        if (width > std::numeric_limits<std::size_t>::max() / height) {
            throw std::overflow_error("Packed texture mip dimensions overflow.");
        }
        const std::size_t texelCount = width * height;
        if (texelCount > maximumUint32 || addedTexels > maximumUint32 - texelCount ||
            scene.textureTexelsRgba8.size() > maximumUint32 - addedTexels - texelCount) {
            throw std::overflow_error("Packed scene texture texels exceed the 32-bit range.");
        }
        const std::size_t channels = static_cast<std::size_t>(image.channels());
        if (texelCount > std::numeric_limits<std::size_t>::max() / channels ||
            image.mipPixels(level).size() != texelCount * channels) {
            throw std::logic_error("A material texture mip has an inconsistent byte count.");
        }
        addedTexels += texelCount;
    }

    const std::size_t newMipCount = scene.textureMipLevels.size() + levelCount;
    const std::size_t newTexelCount = scene.textureTexelsRgba8.size() + addedTexels;
    if (newMipCount > scene.textureMipLevels.max_size() ||
        newTexelCount > scene.textureTexelsRgba8.max_size() ||
        newTexelCount > std::numeric_limits<std::size_t>::max() / sizeof(std::uint32_t)) {
        throw std::overflow_error("Packed texture storage exceeds the addressable range.");
    }

    const std::uint32_t textureId = static_cast<std::uint32_t>(scene.textures.size());
    scene.textures.push_back(PackedTexture{
        static_cast<std::uint32_t>(scene.textureMipLevels.size()),
        static_cast<std::uint32_t>(levelCount),
        static_cast<std::uint32_t>(image.width()),
        static_cast<std::uint32_t>(image.height()),
    });

    for (std::size_t level = 0U; level < levelCount; ++level) {
        const std::uint32_t width = static_cast<std::uint32_t>(image.mipWidth(level));
        const std::uint32_t height = static_cast<std::uint32_t>(image.mipHeight(level));
        const std::size_t texelCount =
            static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
        scene.textureMipLevels.push_back(PackedTextureMip{
            static_cast<std::uint32_t>(scene.textureTexelsRgba8.size()),
            static_cast<std::uint32_t>(texelCount),
            width,
            height,
        });
        const std::vector<std::uint8_t>& pixels = image.mipPixels(level);
        for (std::size_t texel = 0U; texel < texelCount; ++texel) {
            scene.textureTexelsRgba8.push_back(
                packRgba8(pixels, texel, image.channels()));
        }
    }
    return textureId;
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
        triangleCount() > kInvalidSceneId || textures.size() > kInvalidSceneId ||
        textureMipLevels.size() > kInvalidSceneId ||
        textureTexelsRgba8.size() > kInvalidSceneId) {
        throw std::invalid_argument("Packed scene counts exceed the 32-bit index range.");
    }
    if (textureTexelsRgba8.size() >
        std::numeric_limits<std::size_t>::max() / sizeof(std::uint32_t)) {
        throw std::invalid_argument("Packed scene texture byte size overflows.");
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

    std::size_t expectedMipLevel = 0U;
    std::size_t expectedTexelOffset = 0U;
    for (const PackedTexture& texture : textures) {
        if (texture.width == 0U || texture.height == 0U ||
            texture.mipLevelCount == 0U ||
            texture.mipLevelCount > kPackedTextureMaxMipLevels) {
            throw std::invalid_argument("Packed texture metadata is invalid.");
        }
        if (static_cast<std::size_t>(texture.firstMipLevel) != expectedMipLevel ||
            expectedMipLevel > textureMipLevels.size() ||
            static_cast<std::size_t>(texture.mipLevelCount) >
                textureMipLevels.size() - expectedMipLevel) {
            throw std::invalid_argument(
                "Packed texture mip ranges must be contiguous and in range.");
        }

        std::uint32_t expectedWidth = texture.width;
        std::uint32_t expectedHeight = texture.height;
        for (std::uint32_t level = 0U; level < texture.mipLevelCount; ++level) {
            const PackedTextureMip& mip =
                textureMipLevels[expectedMipLevel + static_cast<std::size_t>(level)];
            if (mip.width != expectedWidth || mip.height != expectedHeight) {
                throw std::invalid_argument(
                    "Packed texture mip extents do not form the expected ceil-halved chain.");
            }
            if (static_cast<std::size_t>(mip.texelOffset) != expectedTexelOffset) {
                throw std::invalid_argument(
                    "Packed texture texel ranges must be contiguous.");
            }
            if (level + 1U < texture.mipLevelCount &&
                expectedWidth == 1U && expectedHeight == 1U) {
                throw std::invalid_argument(
                    "Packed textures must not repeat the terminal 1x1 mip.");
            }

            const std::size_t width = static_cast<std::size_t>(mip.width);
            const std::size_t height = static_cast<std::size_t>(mip.height);
            if (width > std::numeric_limits<std::size_t>::max() / height) {
                throw std::invalid_argument("Packed texture mip dimensions overflow.");
            }
            const std::size_t texelCount = width * height;
            if (texelCount > std::numeric_limits<std::uint32_t>::max() ||
                static_cast<std::size_t>(mip.texelCount) != texelCount ||
                expectedTexelOffset > textureTexelsRgba8.size() ||
                texelCount > textureTexelsRgba8.size() - expectedTexelOffset) {
                throw std::invalid_argument("Packed texture mip texel range is invalid.");
            }
            expectedTexelOffset += texelCount;
            expectedWidth = expectedWidth / 2U + expectedWidth % 2U;
            expectedHeight = expectedHeight / 2U + expectedHeight % 2U;
        }

        const PackedTextureMip& lastMip =
            textureMipLevels[
                expectedMipLevel + static_cast<std::size_t>(texture.mipLevelCount) - 1U];
        if (lastMip.width != 1U || lastMip.height != 1U) {
            throw std::invalid_argument("Packed textures must contain a complete mip chain.");
        }
        expectedMipLevel += static_cast<std::size_t>(texture.mipLevelCount);
    }
    if (expectedMipLevel != textureMipLevels.size() ||
        expectedTexelOffset != textureTexelsRgba8.size()) {
        throw std::invalid_argument("Packed texture storage contains unreferenced data.");
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
        for (std::size_t slot = 0U; slot < kTextureFlags.size(); ++slot) {
            const std::uint32_t textureId = material.textureIds[slot];
            const bool hasTextureFlag =
                (material.flagsAndReserved[0] & kTextureFlags[slot]) != 0U;
            const bool hasTextureId = textureId != kInvalidSceneId;
            if (hasTextureFlag != hasTextureId) {
                throw std::invalid_argument(
                    "Packed material texture flags and IDs are inconsistent.");
            }
            if (hasTextureId && textureId >= textures.size()) {
                throw std::invalid_argument(
                    "Packed material contains an out-of-range texture ID.");
            }
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
    std::map<std::string, std::uint32_t> textureIdsByPath;
    for (const Material& material : materials_) {
        PackedMaterial packedMaterial = packMaterial(material);
        for (std::size_t slotIndex = 0U; slotIndex < kTextureSlots.size(); ++slotIndex) {
            const TextureSlot slot = kTextureSlots[slotIndex];
            if (!material.hasTexture(slot)) {
                continue;
            }
            const std::string key =
                normalizedTextureKey(material.textureSourcePath(slot));
            const auto existing = textureIdsByPath.find(key);
            std::uint32_t textureId = kInvalidSceneId;
            if (existing != textureIdsByPath.end()) {
                textureId = existing->second;
            } else {
                textureId = appendTexture(result, material.textureData(slot));
                textureIdsByPath.emplace(key, textureId);
            }
            packedMaterial.textureIds[slotIndex] = textureId;
        }
        result.materials.push_back(packedMaterial);
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
