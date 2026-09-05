#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "raym0nade/geometry.hpp"

namespace raym0nade {

inline constexpr std::size_t kMaxMipmapLevels = 32;

enum class TextureSlot : std::size_t {
    diffuse,
    specular,
    emissive,
    normal,
    count,
};

class ImageData {
public:
    void setBaseLevel(int imageWidth, int imageHeight, int channelCount, std::vector<std::uint8_t> pixels);
    void generateMipmaps();

    [[nodiscard]] vec3 sampleRgb(float u, float v, float depth) const;
    [[nodiscard]] vec4 sampleRgba(float u, float v, float depth) const;

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] bool hasTransparency() const noexcept;
    [[nodiscard]] int width() const noexcept;
    [[nodiscard]] int height() const noexcept;
    [[nodiscard]] int channels() const noexcept;
    [[nodiscard]] std::size_t mipLevelCount() const noexcept;
    [[nodiscard]] int mipWidth(std::size_t level) const;
    [[nodiscard]] int mipHeight(std::size_t level) const;
    [[nodiscard]] const std::vector<std::uint8_t>& mipPixels(std::size_t level) const;

private:
    int width_{0};
    int height_{0};
    int channels_{0};
    std::array<std::vector<std::uint8_t>, kMaxMipmapLevels> levels_;
    std::size_t mipLevelCount_{0};
};

class Material {
public:
    std::string name;
    int id{0};
    bool hasCutoutTransparency{false};
    float opacity{1.0F};
    float ior{1.0F};
    float roughness{0.8F};
    float metallic{0.0F};
    vec3 diffuseFactor{1.0F};
    vec3 emissiveFactor{0.0F};
    vec3 transmissionColor{1.0F};

    void loadTexture(TextureSlot slot, const std::filesystem::path& filename);

    [[nodiscard]] bool hasTexture(TextureSlot slot) const noexcept;
    [[nodiscard]] const ImageData& textureData(TextureSlot slot) const;
    [[nodiscard]] const std::filesystem::path& textureSourcePath(TextureSlot slot) const;
    [[nodiscard]] bool isEmissive() const noexcept;
    [[nodiscard]] vec4 diffuseColor(float u, float v, float footprint) const;
    [[nodiscard]] vec3 normal(float u, float v, float footprint) const;
    [[nodiscard]] vec3 emissiveColor(float u, float v, float footprint) const;
    void surfaceParameters(
        float u, float v, float& surfaceRoughness, float& surfaceMetallic) const;

private:
    static ImageData loadDds(const std::filesystem::path& filename);
    static ImageData loadPng(const std::filesystem::path& filename);
    [[nodiscard]] const ImageData& texture(TextureSlot slot) const;
    [[nodiscard]] ImageData& texture(TextureSlot slot);

    std::array<ImageData, static_cast<std::size_t>(TextureSlot::count)> textures_;
    std::array<std::filesystem::path, static_cast<std::size_t>(TextureSlot::count)>
        textureSourcePaths_;
};

[[nodiscard]] std::string urlDecode(const std::string& source);

}  // namespace raym0nade
