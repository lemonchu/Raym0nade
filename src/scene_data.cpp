#include "raym0nade/scene_data.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
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
    if (values.size() >
        std::numeric_limits<std::uintptr_t>::max() / sizeof(Value)) {
        throw std::logic_error(
            std::string{"Packed scene contains an oversized "} + description + '.');
    }
    const std::uintptr_t beginAddress =
        reinterpret_cast<std::uintptr_t>(values.data());
    const std::uintptr_t pointerAddress =
        reinterpret_cast<std::uintptr_t>(pointer);
    const std::uintptr_t storageBytes =
        static_cast<std::uintptr_t>(values.size()) * sizeof(Value);
    if (pointerAddress < beginAddress) {
        throw std::logic_error(
            std::string{"Packed scene contains an out-of-range "} + description + '.');
    }
    const std::uintptr_t byteOffset = pointerAddress - beginAddress;
    if (byteOffset >= storageBytes || byteOffset % sizeof(Value) != 0U) {
        throw std::logic_error(
            std::string{"Packed scene contains an out-of-range "} + description + '.');
    }
    const std::size_t index =
        static_cast<std::size_t>(byteOffset / sizeof(Value));
    if (values.data() + index != pointer) {
        throw std::logic_error(
            std::string{"Packed scene contains an invalid "} + description + '.');
    }
    return index;
}

bool allFinite(const std::array<float, 4>& values) noexcept {
    return std::all_of(values.begin(), values.end(), [](float value) {
        return std::isfinite(value);
    });
}

float faceArea(const Face& face) noexcept {
    return 0.5F * glm::length(glm::cross(
                      face.vertices[1] - face.vertices[0],
                      face.vertices[2] - face.vertices[0]));
}

float environmentWeight(const vec3& radiance, float solidAngle) noexcept {
    const float luminance = glm::dot(radiance, kLuminanceWeights);
    if (!std::isfinite(luminance) || !(luminance > 0.0F)) {
        return 0.0F;
    }
    const float weight = luminance * solidAngle;
    return std::isfinite(weight) && weight > 0.0F ? weight : 0.0F;
}

bool nearlyEqual(float left, float right, float scale = 1.0F) noexcept {
    if (!std::isfinite(left) || !std::isfinite(right) || !std::isfinite(scale) ||
        !(scale > 0.0F)) {
        return false;
    }
    constexpr float kToleranceInUlps = 16.0F;
    const float magnitude = std::max({std::abs(left), std::abs(right), scale});
    return std::abs(left - right) <=
           kToleranceInUlps * std::numeric_limits<float>::epsilon() * magnitude;
}

vec3 packedPosition(const PackedVertex& vertex) noexcept {
    return vec3{
        vertex.positionAndNormalX[0],
        vertex.positionAndNormalX[1],
        vertex.positionAndNormalX[2],
    };
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

void appendAreaLights(
    PackedSceneData& scene,
    const std::vector<LightObject>& lights,
    const std::vector<VertexData>& sourceVertices,
    const std::vector<Material>& sourceMaterials,
    const std::vector<std::uint32_t>& sourceToPacked) {
    if (lights.size() >= kInvalidSceneId) {
        throw std::overflow_error("Packed scene has too many area lights.");
    }
    scene.areaLights.reserve(lights.size());

    for (const LightObject& light : lights) {
        if (!isFinite(light.center) || !std::isfinite(light.power) ||
            !(light.power > 0.0F) || light.faces.empty() ||
            light.faceDistribution.size() != light.faces.size()) {
            throw std::logic_error("A model area light is structurally invalid.");
        }
        if (light.faces.size() >= kInvalidSceneId ||
            scene.areaLightTriangles.size() >
                static_cast<std::size_t>(kInvalidSceneId) - light.faces.size()) {
            throw std::overflow_error(
                "Packed scene area-light triangles exceed the 32-bit range.");
        }

        const std::uint32_t firstTriangle =
            static_cast<std::uint32_t>(scene.areaLightTriangles.size());
        double cumulativeProbability = 0.0;
        for (std::size_t faceIndex = 0U; faceIndex < light.faces.size(); ++faceIndex) {
            const Face& face = light.faces[faceIndex];
            const float area = faceArea(face);
            const float probability = light.faceDistribution.pdf(faceIndex);
            if (!std::isfinite(area) || !(area > 0.0F) ||
                !std::isfinite(probability) || !(probability > 0.0F)) {
                throw std::logic_error(
                    "A model area-light triangle has an invalid area or probability.");
            }

            PackedAreaLightTriangle packedTriangle;
            for (std::size_t corner = 0U; corner < 3U; ++corner) {
                const std::size_t sourceIndex = borrowedIndex(
                    face.vertexData[corner], sourceVertices, "area-light vertex pointer");
                const std::uint32_t packedIndex = sourceToPacked[sourceIndex];
                if (packedIndex == kInvalidSceneId || packedIndex >= scene.vertices.size() ||
                    !samePosition(scene.vertices[packedIndex], face.vertices[corner])) {
                    throw std::logic_error(
                        "An area-light vertex does not map to packed geometry.");
                }
                packedTriangle.vertexIdsAndMaterialId[corner] = packedIndex;
            }
            packedTriangle.vertexIdsAndMaterialId[3] = static_cast<std::uint32_t>(
                borrowedIndex(face.material, sourceMaterials, "area-light material pointer"));

            cumulativeProbability += static_cast<double>(probability);
            const float cdf =
                faceIndex + 1U == light.faces.size()
                    ? 1.0F
                    : std::clamp(
                          static_cast<float>(cumulativeProbability), 0.0F, 1.0F);
            packedTriangle.areaProbabilityCdfAndReserved = {
                area,
                probability,
                cdf,
                0.0F,
            };
            scene.areaLightTriangles.push_back(packedTriangle);
        }

        PackedAreaLight packedLight;
        packedLight.centerAndPower = {
            light.center.x,
            light.center.y,
            light.center.z,
            light.power,
        };
        packedLight.triangleRangeAndReserved = {
            firstTriangle,
            static_cast<std::uint32_t>(light.faces.size()),
            0U,
            0U,
        };
        scene.areaLights.push_back(packedLight);
    }
}

void appendEnvironment(PackedSceneData& scene, const SkyBox& sky) {
    if (sky.empty()) {
        return;
    }
    if (sky.width() <= 0 || sky.height() <= 0) {
        throw std::logic_error("A non-empty environment has an invalid extent.");
    }

    const std::size_t width = static_cast<std::size_t>(sky.width());
    const std::size_t height = static_cast<std::size_t>(sky.height());
    if (width > std::numeric_limits<std::size_t>::max() / height) {
        throw std::overflow_error("Packed environment dimensions overflow.");
    }
    const std::size_t texelCount = width * height;
    const std::vector<vec3>& radiance = sky.radiancePixels();
    if (texelCount >= kInvalidSceneId || radiance.size() != texelCount) {
        throw std::logic_error(
            "Packed environment radiance does not match its 32-bit extent.");
    }

    std::vector<float> solidAngles(height, 0.0F);
    std::vector<double> rowWeights(height, 0.0);
    double totalWeight = 0.0;
    const double azimuthWidth =
        2.0 * static_cast<double>(kPi) / static_cast<double>(width);
    for (std::size_t row = 0U; row < height; ++row) {
        const double polarMinimum =
            static_cast<double>(kPi) * static_cast<double>(row) /
            static_cast<double>(height);
        const double polarMaximum =
            static_cast<double>(kPi) * static_cast<double>(row + 1U) /
            static_cast<double>(height);
        const float solidAngle = static_cast<float>(
            azimuthWidth * (std::cos(polarMinimum) - std::cos(polarMaximum)));
        if (!std::isfinite(solidAngle) || !(solidAngle > 0.0F)) {
            throw std::logic_error("Packed environment row has an invalid solid angle.");
        }
        solidAngles[row] = solidAngle;

        for (std::size_t column = 0U; column < width; ++column) {
            const vec3& value = radiance[row * width + column];
            if (!isFinite(value) || glm::any(glm::lessThan(value, vec3{0.0F}))) {
                throw std::logic_error(
                    "Packed environment radiance must be finite and non-negative.");
            }
            const double weight =
                static_cast<double>(environmentWeight(value, solidAngle));
            rowWeights[row] += weight;
            totalWeight += weight;
        }
    }

    const bool hasImportance =
        std::isfinite(totalWeight) && totalWeight > 0.0;
    scene.environment = PackedEnvironment{
        static_cast<std::uint32_t>(width),
        static_cast<std::uint32_t>(height),
        hasImportance ? kPackedEnvironmentHasImportance : 0U,
        0U,
    };
    scene.environmentRows.reserve(height);
    scene.environmentTexels.reserve(texelCount);

    double cumulativeRowWeight = 0.0;
    for (std::size_t row = 0U; row < height; ++row) {
        const bool rowHasImportance =
            hasImportance && std::isfinite(rowWeights[row]) && rowWeights[row] > 0.0;
        const float rowProbability = rowHasImportance
                                         ? static_cast<float>(rowWeights[row] / totalWeight)
                                         : 0.0F;
        cumulativeRowWeight += rowHasImportance ? rowWeights[row] : 0.0;
        const float rowCdf =
            !hasImportance
                ? 0.0F
                : (row + 1U == height
                       ? 1.0F
                       : std::clamp(
                             static_cast<float>(cumulativeRowWeight / totalWeight),
                             0.0F,
                             1.0F));
        scene.environmentRows.push_back(PackedEnvironmentRow{
            rowProbability,
            rowCdf,
            solidAngles[row],
            0.0F,
        });

        double cumulativeColumnWeight = 0.0;
        for (std::size_t column = 0U; column < width; ++column) {
            const vec3& value = radiance[row * width + column];
            cumulativeColumnWeight +=
                rowHasImportance
                    ? static_cast<double>(environmentWeight(value, solidAngles[row]))
                    : 0.0;
            const float conditionalCdf =
                column + 1U == width
                    ? 1.0F
                    : (rowHasImportance
                           ? std::clamp(
                                 static_cast<float>(
                                     cumulativeColumnWeight / rowWeights[row]),
                                 0.0F,
                                 1.0F)
                           : static_cast<float>(
                                 static_cast<double>(column + 1U) /
                                 static_cast<double>(width)));
            scene.environmentTexels.push_back(PackedEnvironmentTexel{{
                value.x,
                value.y,
                value.z,
                conditionalCdf,
            }});
        }
    }
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
        textureTexelsRgba8.size() > kInvalidSceneId ||
        areaLights.size() > kInvalidSceneId ||
        areaLightTriangles.size() > kInvalidSceneId ||
        environmentRows.size() > kInvalidSceneId ||
        environmentTexels.size() > kInvalidSceneId) {
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

    std::size_t expectedAreaTriangle = 0U;
    for (const PackedAreaLight& light : areaLights) {
        if (!allFinite(light.centerAndPower) || !(light.centerAndPower[3] > 0.0F) ||
            light.triangleRangeAndReserved[1] == 0U ||
            light.triangleRangeAndReserved[2] != 0U ||
            light.triangleRangeAndReserved[3] != 0U) {
            throw std::invalid_argument("Packed area-light metadata is invalid.");
        }
        if (static_cast<std::size_t>(light.triangleRangeAndReserved[0]) !=
                expectedAreaTriangle ||
            expectedAreaTriangle > areaLightTriangles.size() ||
            static_cast<std::size_t>(light.triangleRangeAndReserved[1]) >
                areaLightTriangles.size() - expectedAreaTriangle) {
            throw std::invalid_argument(
                "Packed area-light triangle ranges must be contiguous and in range.");
        }

        float previousCdf = 0.0F;
        const std::size_t triangleCount =
            static_cast<std::size_t>(light.triangleRangeAndReserved[1]);
        for (std::size_t localIndex = 0U; localIndex < triangleCount; ++localIndex) {
            const PackedAreaLightTriangle& triangle =
                areaLightTriangles[expectedAreaTriangle + localIndex];
            for (std::size_t corner = 0U; corner < 3U; ++corner) {
                if (triangle.vertexIdsAndMaterialId[corner] >= vertices.size()) {
                    throw std::invalid_argument(
                        "Packed area light contains an out-of-range vertex ID.");
                }
            }
            const std::uint32_t materialId = triangle.vertexIdsAndMaterialId[3];
            if (materialId >= materials.size()) {
                throw std::invalid_argument(
                    "Packed area light contains an out-of-range material ID.");
            }
            const PackedMaterial& material = materials[materialId];
            const bool hasEmission =
                material.emissionAndIor[0] > 0.0F ||
                material.emissionAndIor[1] > 0.0F ||
                material.emissionAndIor[2] > 0.0F ||
                (material.flagsAndReserved[0] &
                 kPackedMaterialHasEmissiveTexture) != 0U;
            if (!hasEmission) {
                throw std::invalid_argument(
                    "Packed area-light triangles must reference emissive materials.");
            }

            const std::array<float, 4>& values =
                triangle.areaProbabilityCdfAndReserved;
            if (!allFinite(values) || !(values[0] > 0.0F) ||
                !(values[1] > 0.0F) || values[1] > 1.0F ||
                values[2] < previousCdf || values[2] > 1.0F ||
                values[3] != 0.0F) {
                throw std::invalid_argument(
                    "Packed area-light triangle sampling data is invalid.");
            }
            const vec3 position0 =
                packedPosition(vertices[triangle.vertexIdsAndMaterialId[0]]);
            const vec3 position1 =
                packedPosition(vertices[triangle.vertexIdsAndMaterialId[1]]);
            const vec3 position2 =
                packedPosition(vertices[triangle.vertexIdsAndMaterialId[2]]);
            const float geometricArea = 0.5F * glm::length(
                glm::cross(position1 - position0, position2 - position0));
            if (!(geometricArea > 0.0F) || !std::isfinite(geometricArea) ||
                !nearlyEqual(
                    values[0],
                    geometricArea,
                    std::numeric_limits<float>::min())) {
                throw std::invalid_argument(
                    "Packed area-light area does not match its triangle.");
            }
            const float cdfProbability = values[2] - previousCdf;
            if (!nearlyEqual(cdfProbability, values[1])) {
                throw std::invalid_argument(
                    "Packed area-light face probability and CDF are inconsistent.");
            }
            previousCdf = values[2];
        }
        if (previousCdf != 1.0F) {
            throw std::invalid_argument(
                "Packed area-light face CDF must end at one.");
        }
        expectedAreaTriangle += triangleCount;
    }
    if (expectedAreaTriangle != areaLightTriangles.size()) {
        throw std::invalid_argument(
            "Packed area-light storage contains unreferenced triangles.");
    }

    if ((environment.flags & ~kPackedEnvironmentKnownFlags) != 0U ||
        environment.reserved != 0U) {
        throw std::invalid_argument("Packed environment metadata is invalid.");
    }
    const bool hasEnvironmentExtent =
        environment.width != 0U || environment.height != 0U;
    if (!hasEnvironmentExtent) {
        if (environment.width != 0U || environment.height != 0U ||
            environment.flags != 0U || !environmentRows.empty() ||
            !environmentTexels.empty()) {
            throw std::invalid_argument(
                "An empty packed environment must not own sampling data.");
        }
        return;
    }
    if (environment.width == 0U || environment.height == 0U) {
        throw std::invalid_argument(
            "Packed environment width and height must both be nonzero.");
    }
    const std::size_t environmentWidth =
        static_cast<std::size_t>(environment.width);
    const std::size_t environmentHeight =
        static_cast<std::size_t>(environment.height);
    if (environmentWidth >
        std::numeric_limits<std::size_t>::max() / environmentHeight) {
        throw std::invalid_argument("Packed environment dimensions overflow.");
    }
    const std::size_t environmentTexelCount =
        environmentWidth * environmentHeight;
    if (environmentRows.size() != environmentHeight ||
        environmentTexels.size() != environmentTexelCount) {
        throw std::invalid_argument(
            "Packed environment arrays do not match its extent.");
    }

    std::vector<float> expectedSolidAngles(environmentHeight, 0.0F);
    std::vector<double> expectedRowWeights(environmentHeight, 0.0);
    double expectedTotalWeight = 0.0;
    const double azimuthWidth =
        2.0 * static_cast<double>(kPi) /
        static_cast<double>(environmentWidth);
    for (std::size_t rowIndex = 0U; rowIndex < environmentHeight; ++rowIndex) {
        const double polarMinimum =
            static_cast<double>(kPi) * static_cast<double>(rowIndex) /
            static_cast<double>(environmentHeight);
        const double polarMaximum =
            static_cast<double>(kPi) * static_cast<double>(rowIndex + 1U) /
            static_cast<double>(environmentHeight);
        const float expectedSolidAngle = static_cast<float>(
            azimuthWidth *
            (std::cos(polarMinimum) - std::cos(polarMaximum)));
        if (!std::isfinite(expectedSolidAngle) ||
            !(expectedSolidAngle > 0.0F)) {
            throw std::invalid_argument(
                "Packed environment dimensions produce an invalid solid angle.");
        }
        expectedSolidAngles[rowIndex] = expectedSolidAngle;
        for (std::size_t column = 0U; column < environmentWidth; ++column) {
            const PackedEnvironmentTexel& texel =
                environmentTexels[rowIndex * environmentWidth + column];
            if (!allFinite(texel.radianceAndConditionalCdf) ||
                texel.radianceAndConditionalCdf[0] < 0.0F ||
                texel.radianceAndConditionalCdf[1] < 0.0F ||
                texel.radianceAndConditionalCdf[2] < 0.0F) {
                throw std::invalid_argument(
                    "Packed environment radiance must be finite and non-negative.");
            }
            const vec3 radiance{
                texel.radianceAndConditionalCdf[0],
                texel.radianceAndConditionalCdf[1],
                texel.radianceAndConditionalCdf[2],
            };
            const double weight =
                static_cast<double>(
                    environmentWeight(radiance, expectedSolidAngle));
            expectedRowWeights[rowIndex] += weight;
            expectedTotalWeight += weight;
        }
    }

    const bool hasImportance =
        (environment.flags & kPackedEnvironmentHasImportance) != 0U;
    const bool expectsImportance =
        std::isfinite(expectedTotalWeight) && expectedTotalWeight > 0.0;
    if (hasImportance != expectsImportance) {
        throw std::invalid_argument(
            "Packed environment importance metadata does not match its radiance.");
    }

    double cumulativeRowWeight = 0.0;
    for (std::size_t rowIndex = 0U; rowIndex < environmentHeight; ++rowIndex) {
        const PackedEnvironmentRow& row = environmentRows[rowIndex];
        if (!std::isfinite(row.probability) ||
            !std::isfinite(row.cumulativeProbability) ||
            !std::isfinite(row.solidAngle) || row.probability < 0.0F ||
            row.probability > 1.0F || row.cumulativeProbability < 0.0F ||
            row.cumulativeProbability > 1.0F ||
            !(row.solidAngle > 0.0F) ||
            row.reserved != 0.0F) {
            throw std::invalid_argument(
                "Packed environment row metadata is invalid.");
        }
        if (!nearlyEqual(
                row.solidAngle,
                expectedSolidAngles[rowIndex],
                std::numeric_limits<float>::min())) {
            throw std::invalid_argument(
                "Packed environment row solid angle does not match its extent.");
        }

        const bool rowHasImportance =
            hasImportance && std::isfinite(expectedRowWeights[rowIndex]) &&
            expectedRowWeights[rowIndex] > 0.0;
        const float expectedRowProbability =
            rowHasImportance
                ? static_cast<float>(
                      expectedRowWeights[rowIndex] / expectedTotalWeight)
                : 0.0F;
        cumulativeRowWeight +=
            rowHasImportance ? expectedRowWeights[rowIndex] : 0.0;
        const float expectedRowCdf =
            !hasImportance
                ? 0.0F
                : (rowIndex + 1U == environmentHeight
                       ? 1.0F
                       : std::clamp(
                             static_cast<float>(
                                 cumulativeRowWeight / expectedTotalWeight),
                             0.0F,
                             1.0F));
        if (!nearlyEqual(row.probability, expectedRowProbability) ||
            !nearlyEqual(row.cumulativeProbability, expectedRowCdf)) {
            throw std::invalid_argument(
                "Packed environment row distribution does not match its radiance.");
        }

        double cumulativeColumnWeight = 0.0;
        for (std::size_t column = 0U; column < environmentWidth; ++column) {
            const PackedEnvironmentTexel& texel =
                environmentTexels[rowIndex * environmentWidth + column];
            if (texel.radianceAndConditionalCdf[3] < 0.0F ||
                texel.radianceAndConditionalCdf[3] > 1.0F) {
                throw std::invalid_argument(
                    "Packed environment conditional CDF is outside the unit interval.");
            }
            const vec3 radiance{
                texel.radianceAndConditionalCdf[0],
                texel.radianceAndConditionalCdf[1],
                texel.radianceAndConditionalCdf[2],
            };
            cumulativeColumnWeight +=
                rowHasImportance
                    ? static_cast<double>(
                          environmentWeight(
                              radiance, expectedSolidAngles[rowIndex]))
                    : 0.0;
            const float expectedConditionalCdf =
                column + 1U == environmentWidth
                    ? 1.0F
                    : (rowHasImportance
                           ? std::clamp(
                                 static_cast<float>(
                                     cumulativeColumnWeight /
                                     expectedRowWeights[rowIndex]),
                                 0.0F,
                                 1.0F)
                           : static_cast<float>(
                                 static_cast<double>(column + 1U) /
                                 static_cast<double>(environmentWidth)));
            if (!nearlyEqual(
                    texel.radianceAndConditionalCdf[3],
                    expectedConditionalCdf)) {
                throw std::invalid_argument(
                    "Packed environment conditional CDF does not match its radiance.");
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

    appendAreaLights(result, lights_, vertexData_, materials_, sourceToPacked);
    appendEnvironment(result, sky_);
    result.validate();
    return result;
}

}  // namespace raym0nade
