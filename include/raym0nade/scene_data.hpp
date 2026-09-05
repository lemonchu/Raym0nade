#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace raym0nade {

inline constexpr std::uint32_t kPackedSceneFormatVersion = 4U;
inline constexpr std::uint32_t kInvalidSceneId = 0xffffffffU;
inline constexpr std::uint32_t kPackedTextureMaxMipLevels = 32U;
inline constexpr std::uint32_t kPackedEnvironmentHasImportance = 1U << 0U;
inline constexpr std::uint32_t kPackedEnvironmentKnownFlags =
    kPackedEnvironmentHasImportance;
inline constexpr std::uint32_t kPackedMaterialCutout = 1U << 0U;
inline constexpr std::uint32_t kPackedMaterialHasDiffuseTexture = 1U << 1U;
inline constexpr std::uint32_t kPackedMaterialHasSpecularTexture = 1U << 2U;
inline constexpr std::uint32_t kPackedMaterialHasEmissiveTexture = 1U << 3U;
inline constexpr std::uint32_t kPackedMaterialHasNormalTexture = 1U << 4U;
inline constexpr std::uint32_t kPackedMaterialKnownFlags =
    kPackedMaterialCutout | kPackedMaterialHasDiffuseTexture |
    kPackedMaterialHasSpecularTexture | kPackedMaterialHasEmissiveTexture |
    kPackedMaterialHasNormalTexture;

struct alignas(16) PackedVertex {
    std::array<float, 4> positionAndNormalX{};
    std::array<float, 4> normalYZAndUv{};
};

struct alignas(16) PackedMaterial {
    std::array<float, 4> diffuseAndOpacity{};
    std::array<float, 4> emissionAndIor{};
    std::array<float, 4> transmissionAndRoughness{};
    std::array<float, 4> metallicSpecularAndReserved{};
    std::array<std::uint32_t, 4> textureIds{
        kInvalidSceneId,
        kInvalidSceneId,
        kInvalidSceneId,
        kInvalidSceneId,
    };
    std::array<std::uint32_t, 4> flagsAndReserved{};
};

struct alignas(16) PackedTexture {
    std::uint32_t firstMipLevel{0U};
    std::uint32_t mipLevelCount{0U};
    std::uint32_t width{0U};
    std::uint32_t height{0U};
};

struct alignas(16) PackedTextureMip {
    std::uint32_t texelOffset{0U};
    std::uint32_t texelCount{0U};
    std::uint32_t width{0U};
    std::uint32_t height{0U};
};

struct alignas(16) PackedAreaLight {
    std::array<float, 4> centerAndPower{};
    // first triangle, triangle count, reserved, reserved
    std::array<std::uint32_t, 4> triangleRangeAndReserved{};
};

struct alignas(16) PackedAreaLightTriangle {
    // Three packed vertex IDs followed by the material ID.
    std::array<std::uint32_t, 4> vertexIdsAndMaterialId{};
    // Area, normalized face probability, normalized face CDF, reserved.
    std::array<float, 4> areaProbabilityCdfAndReserved{};
};

struct alignas(16) PackedEnvironment {
    std::uint32_t width{0U};
    std::uint32_t height{0U};
    std::uint32_t flags{0U};
    std::uint32_t reserved{0U};
};

struct alignas(16) PackedEnvironmentRow {
    float probability{0.0F};
    float cumulativeProbability{0.0F};
    float solidAngle{0.0F};
    float reserved{0.0F};
};

struct alignas(16) PackedEnvironmentTexel {
    // Linear RGB radiance followed by the normalized conditional CDF within its row.
    std::array<float, 4> radianceAndConditionalCdf{};
};

class PackedSceneData {
public:
    std::uint32_t formatVersion{kPackedSceneFormatVersion};
    std::vector<PackedVertex> vertices;
    std::vector<std::uint32_t> triangleIndices;
    std::vector<std::uint32_t> triangleMaterialIds;
    std::vector<PackedMaterial> materials;
    std::vector<PackedTexture> textures;
    std::vector<PackedTextureMip> textureMipLevels;
    // Encoded UNORM color/data values in a linear address layout. Each word stores
    // R, G, B, and A in bits 0-7, 8-15, 16-23, and 24-31 respectively.
    std::vector<std::uint32_t> textureTexelsRgba8;
    std::vector<PackedAreaLight> areaLights;
    std::vector<PackedAreaLightTriangle> areaLightTriangles;
    PackedEnvironment environment;
    std::vector<PackedEnvironmentRow> environmentRows;
    std::vector<PackedEnvironmentTexel> environmentTexels;

    [[nodiscard]] std::size_t triangleCount() const noexcept;
    void validate() const;
};

static_assert(sizeof(PackedVertex) == 32U, "PackedVertex must match the GPU ABI.");
static_assert(alignof(PackedVertex) == 16U, "PackedVertex must remain 16-byte aligned.");
static_assert(offsetof(PackedVertex, positionAndNormalX) == 0U);
static_assert(offsetof(PackedVertex, normalYZAndUv) == 16U);
static_assert(std::is_standard_layout_v<PackedVertex>);
static_assert(std::is_trivially_copyable_v<PackedVertex>);
static_assert(sizeof(PackedMaterial) == 96U, "PackedMaterial must match the GPU ABI.");
static_assert(alignof(PackedMaterial) == 16U, "PackedMaterial must remain 16-byte aligned.");
static_assert(offsetof(PackedMaterial, diffuseAndOpacity) == 0U);
static_assert(offsetof(PackedMaterial, emissionAndIor) == 16U);
static_assert(offsetof(PackedMaterial, transmissionAndRoughness) == 32U);
static_assert(offsetof(PackedMaterial, metallicSpecularAndReserved) == 48U);
static_assert(offsetof(PackedMaterial, textureIds) == 64U);
static_assert(offsetof(PackedMaterial, flagsAndReserved) == 80U);
static_assert(std::is_standard_layout_v<PackedMaterial>);
static_assert(std::is_trivially_copyable_v<PackedMaterial>);
static_assert(sizeof(PackedTexture) == 16U, "PackedTexture must match the GPU ABI.");
static_assert(alignof(PackedTexture) == 16U, "PackedTexture must remain 16-byte aligned.");
static_assert(offsetof(PackedTexture, firstMipLevel) == 0U);
static_assert(offsetof(PackedTexture, mipLevelCount) == 4U);
static_assert(offsetof(PackedTexture, width) == 8U);
static_assert(offsetof(PackedTexture, height) == 12U);
static_assert(std::is_standard_layout_v<PackedTexture>);
static_assert(std::is_trivially_copyable_v<PackedTexture>);
static_assert(sizeof(PackedTextureMip) == 16U, "PackedTextureMip must match the GPU ABI.");
static_assert(alignof(PackedTextureMip) == 16U, "PackedTextureMip must remain 16-byte aligned.");
static_assert(offsetof(PackedTextureMip, texelOffset) == 0U);
static_assert(offsetof(PackedTextureMip, texelCount) == 4U);
static_assert(offsetof(PackedTextureMip, width) == 8U);
static_assert(offsetof(PackedTextureMip, height) == 12U);
static_assert(std::is_standard_layout_v<PackedTextureMip>);
static_assert(std::is_trivially_copyable_v<PackedTextureMip>);
static_assert(sizeof(PackedAreaLight) == 32U, "PackedAreaLight must match the GPU ABI.");
static_assert(alignof(PackedAreaLight) == 16U, "PackedAreaLight must remain 16-byte aligned.");
static_assert(offsetof(PackedAreaLight, centerAndPower) == 0U);
static_assert(offsetof(PackedAreaLight, triangleRangeAndReserved) == 16U);
static_assert(std::is_standard_layout_v<PackedAreaLight>);
static_assert(std::is_trivially_copyable_v<PackedAreaLight>);
static_assert(
    sizeof(PackedAreaLightTriangle) == 32U,
    "PackedAreaLightTriangle must match the GPU ABI.");
static_assert(
    alignof(PackedAreaLightTriangle) == 16U,
    "PackedAreaLightTriangle must remain 16-byte aligned.");
static_assert(offsetof(PackedAreaLightTriangle, vertexIdsAndMaterialId) == 0U);
static_assert(offsetof(PackedAreaLightTriangle, areaProbabilityCdfAndReserved) == 16U);
static_assert(std::is_standard_layout_v<PackedAreaLightTriangle>);
static_assert(std::is_trivially_copyable_v<PackedAreaLightTriangle>);
static_assert(sizeof(PackedEnvironment) == 16U, "PackedEnvironment must match the GPU ABI.");
static_assert(
    alignof(PackedEnvironment) == 16U,
    "PackedEnvironment must remain 16-byte aligned.");
static_assert(offsetof(PackedEnvironment, width) == 0U);
static_assert(offsetof(PackedEnvironment, height) == 4U);
static_assert(offsetof(PackedEnvironment, flags) == 8U);
static_assert(offsetof(PackedEnvironment, reserved) == 12U);
static_assert(std::is_standard_layout_v<PackedEnvironment>);
static_assert(std::is_trivially_copyable_v<PackedEnvironment>);
static_assert(
    sizeof(PackedEnvironmentRow) == 16U,
    "PackedEnvironmentRow must match the GPU ABI.");
static_assert(
    alignof(PackedEnvironmentRow) == 16U,
    "PackedEnvironmentRow must remain 16-byte aligned.");
static_assert(offsetof(PackedEnvironmentRow, probability) == 0U);
static_assert(offsetof(PackedEnvironmentRow, cumulativeProbability) == 4U);
static_assert(offsetof(PackedEnvironmentRow, solidAngle) == 8U);
static_assert(offsetof(PackedEnvironmentRow, reserved) == 12U);
static_assert(std::is_standard_layout_v<PackedEnvironmentRow>);
static_assert(std::is_trivially_copyable_v<PackedEnvironmentRow>);
static_assert(
    sizeof(PackedEnvironmentTexel) == 16U,
    "PackedEnvironmentTexel must match the GPU ABI.");
static_assert(
    alignof(PackedEnvironmentTexel) == 16U,
    "PackedEnvironmentTexel must remain 16-byte aligned.");
static_assert(offsetof(PackedEnvironmentTexel, radianceAndConditionalCdf) == 0U);
static_assert(std::is_standard_layout_v<PackedEnvironmentTexel>);
static_assert(std::is_trivially_copyable_v<PackedEnvironmentTexel>);

}  // namespace raym0nade
