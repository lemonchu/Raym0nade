#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace raym0nade {

inline constexpr std::uint32_t kPackedSceneFormatVersion = 2U;
inline constexpr std::uint32_t kInvalidSceneId = 0xffffffffU;
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

class PackedSceneData {
public:
    std::uint32_t formatVersion{kPackedSceneFormatVersion};
    std::vector<PackedVertex> vertices;
    std::vector<std::uint32_t> triangleIndices;
    std::vector<std::uint32_t> triangleMaterialIds;
    std::vector<PackedMaterial> materials;

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

}  // namespace raym0nade
