#include "raym0nade/material.hpp"

#include <png.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "python_image_loader.hpp"

namespace raym0nade {
namespace {

class FileHandle {
public:
    explicit FileHandle(const std::filesystem::path& filename) {
#if defined(_WIN32) && defined(_MSC_VER)
        if (_wfopen_s(&file_, filename.c_str(), L"rb") != 0) {
            file_ = nullptr;
        }
#elif defined(_WIN32)
        file_ = _wfopen(filename.c_str(), L"rb");
#else
        file_ = std::fopen(filename.c_str(), "rb");
#endif
        if (file_ == nullptr) {
            throw std::runtime_error("Failed to open image: " + filename.string());
        }
    }

    FileHandle(const FileHandle&) = delete;
    FileHandle& operator=(const FileHandle&) = delete;

    ~FileHandle() {
        if (file_ != nullptr) {
            std::fclose(file_);
        }
    }

    [[nodiscard]] FILE* get() const noexcept {
        return file_;
    }

private:
    FILE* file_{nullptr};
};

class PngImageHandle {
public:
    PngImageHandle() {
        image_.version = PNG_IMAGE_VERSION;
    }

    PngImageHandle(const PngImageHandle&) = delete;
    PngImageHandle& operator=(const PngImageHandle&) = delete;

    ~PngImageHandle() {
        png_image_free(&image_);
    }

    [[nodiscard]] png_image& get() noexcept {
        return image_;
    }

private:
    png_image image_{};
};

std::uint8_t component(
    const std::vector<std::uint8_t>& pixels, std::size_t pixelIndex, int channels, int requestedChannel) {
    const std::size_t base = pixelIndex * static_cast<std::size_t>(channels);
    if (channels == 1) {
        return requestedChannel == 3 ? 255 : pixels[base];
    }
    if (channels == 2) {
        return requestedChannel == 3 ? pixels[base + 1] : pixels[base];
    }
    if (requestedChannel == 3) {
        return channels >= 4 ? pixels[base + 3] : 255;
    }
    return pixels[base + static_cast<std::size_t>(requestedChannel)];
}

template <typename Vector>
Vector pixel(const std::vector<std::uint8_t>& pixels, std::size_t pixelIndex, int channels);

template <>
vec3 pixel<vec3>(const std::vector<std::uint8_t>& pixels, std::size_t pixelIndex, int channels) {
    return vec3{
               static_cast<float>(component(pixels, pixelIndex, channels, 0)),
               static_cast<float>(component(pixels, pixelIndex, channels, 1)),
               static_cast<float>(component(pixels, pixelIndex, channels, 2))} /
           255.0F;
}

template <>
vec4 pixel<vec4>(const std::vector<std::uint8_t>& pixels, std::size_t pixelIndex, int channels) {
    return vec4{
               static_cast<float>(component(pixels, pixelIndex, channels, 0)),
               static_cast<float>(component(pixels, pixelIndex, channels, 1)),
               static_cast<float>(component(pixels, pixelIndex, channels, 2)),
               static_cast<float>(component(pixels, pixelIndex, channels, 3))} /
           255.0F;
}

int wrapIndex(int value, int extent) noexcept {
    value %= extent;
    return value < 0 ? value + extent : value;
}

float wrapUnit(float value) noexcept {
    if (!std::isfinite(value)) {
        return 0.0F;
    }
    const float wrapped = value - std::floor(value);
    return wrapped >= 0.0F && wrapped < 1.0F ? wrapped : 0.0F;
}

int mipExtent(int baseExtent, std::size_t level) noexcept {
    int extent = std::max(baseExtent, 1);
    for (std::size_t index = 0; index < level && extent > 1; ++index) {
        extent = extent / 2 + extent % 2;
    }
    return extent;
}

template <typename Vector>
Vector bilinearSample(
    const std::vector<std::uint8_t>& pixels,
    int width,
    int height,
    int channels,
    float u,
    float v) {
    const float x = wrapUnit(u) * static_cast<float>(width);
    const float y = wrapUnit(v) * static_cast<float>(height);
    const int unwrappedX = static_cast<int>(std::floor(x));
    const int unwrappedY = static_cast<int>(std::floor(y));
    const float blendX = x - static_cast<float>(unwrappedX);
    const float blendY = y - static_cast<float>(unwrappedY);

    const int x0 = wrapIndex(unwrappedX, width);
    const int y0 = wrapIndex(unwrappedY, height);
    const int x1 = (x0 + 1) % width;
    const int y1 = (y0 + 1) % height;

    const auto index = [width](int sampleX, int sampleY) {
        return static_cast<std::size_t>(sampleY) * static_cast<std::size_t>(width) +
               static_cast<std::size_t>(sampleX);
    };
    const Vector c00 = pixel<Vector>(pixels, index(x0, y0), channels);
    const Vector c01 = pixel<Vector>(pixels, index(x1, y0), channels);
    const Vector c10 = pixel<Vector>(pixels, index(x0, y1), channels);
    const Vector c11 = pixel<Vector>(pixels, index(x1, y1), channels);
    const Vector top = c00 * (1.0F - blendX) + c01 * blendX;
    const Vector bottom = c10 * (1.0F - blendX) + c11 * blendX;
    return top * (1.0F - blendY) + bottom * blendY;
}

float mipDepth(float footprint, int textureWidth) noexcept {
    if (!isFinite(footprint) || footprint <= 0.0F || textureWidth <= 0) {
        return 0.0F;
    }
    return std::max(0.0F, std::log2(footprint * static_cast<float>(textureWidth)));
}

void convertToLinear(vec4& value) noexcept {
    constexpr float gamma = 2.2F;
    value.r = std::pow(std::max(value.r, 0.0F), gamma);
    value.g = std::pow(std::max(value.g, 0.0F), gamma);
    value.b = std::pow(std::max(value.b, 0.0F), gamma);
}

bool isHexDigit(char value) noexcept {
    return (value >= '0' && value <= '9') || (value >= 'a' && value <= 'f') ||
           (value >= 'A' && value <= 'F');
}

unsigned int hexValue(char value) noexcept {
    if (value >= '0' && value <= '9') {
        return static_cast<unsigned int>(value - '0');
    }
    if (value >= 'a' && value <= 'f') {
        return static_cast<unsigned int>(value - 'a' + 10);
    }
    return static_cast<unsigned int>(value - 'A' + 10);
}

}  // namespace

void ImageData::setBaseLevel(
    int imageWidth, int imageHeight, int channelCount, std::vector<std::uint8_t> pixels) {
    if (imageWidth <= 0 || imageHeight <= 0 || channelCount < 1 || channelCount > 4) {
        throw std::invalid_argument("Invalid image dimensions or channel count.");
    }
    const std::size_t width = static_cast<std::size_t>(imageWidth);
    const std::size_t height = static_cast<std::size_t>(imageHeight);
    const std::size_t channels = static_cast<std::size_t>(channelCount);
    const std::size_t maximumSize = std::numeric_limits<std::size_t>::max();
    if (width > maximumSize / height || width * height > maximumSize / channels) {
        throw std::invalid_argument("Image dimensions overflow the addressable buffer size.");
    }
    const std::size_t expectedSize = width * height * channels;
    if (pixels.size() != expectedSize) {
        throw std::invalid_argument("Image buffer size does not match its dimensions.");
    }

    width_ = imageWidth;
    height_ = imageHeight;
    channels_ = channelCount;
    for (auto& level : levels_) {
        level.clear();
    }
    levels_[0] = std::move(pixels);
    mipLevelCount_ = 1;
}

void ImageData::generateMipmaps() {
    if (empty()) {
        return;
    }

    for (std::size_t level = 1; level < kMaxMipmapLevels; ++level) {
        const int previousWidth = mipExtent(width_, level - 1);
        const int previousHeight = mipExtent(height_, level - 1);
        if (previousWidth == 1 && previousHeight == 1) {
            break;
        }
        const int currentWidth = previousWidth / 2 + previousWidth % 2;
        const int currentHeight = previousHeight / 2 + previousHeight % 2;

        auto& current = levels_[level];
        const auto& previous = levels_[level - 1];
        current.resize(static_cast<std::size_t>(currentWidth) * static_cast<std::size_t>(currentHeight) *
                       static_cast<std::size_t>(channels_));
        for (int y = 0; y < currentHeight; ++y) {
            for (int x = 0; x < currentWidth; ++x) {
                for (int channel = 0; channel < channels_; ++channel) {
                    unsigned int sum = 0;
                    for (int offsetY = 0; offsetY < 2; ++offsetY) {
                        for (int offsetX = 0; offsetX < 2; ++offsetX) {
                            const int sourceX = (x * 2 + offsetX) % previousWidth;
                            const int sourceY = (y * 2 + offsetY) % previousHeight;
                            const std::size_t sourceIndex =
                                (static_cast<std::size_t>(sourceY) * static_cast<std::size_t>(previousWidth) +
                                 static_cast<std::size_t>(sourceX)) *
                                    static_cast<std::size_t>(channels_) +
                                static_cast<std::size_t>(channel);
                            sum += previous[sourceIndex];
                        }
                    }
                    const std::size_t destinationIndex =
                        (static_cast<std::size_t>(y) * static_cast<std::size_t>(currentWidth) +
                         static_cast<std::size_t>(x)) *
                            static_cast<std::size_t>(channels_) +
                        static_cast<std::size_t>(channel);
                    current[destinationIndex] = static_cast<std::uint8_t>(sum / 4U);
                }
            }
        }
        mipLevelCount_ = level + 1;
    }
}

namespace {

template <typename Vector>
Vector sampleImage(
    const std::array<std::vector<std::uint8_t>, kMaxMipmapLevels>& levels,
    std::size_t mipLevelCount,
    int width,
    int height,
    int channels,
    float u,
    float v,
    float depth) {
    if (mipLevelCount == 0 || levels[0].empty() || !isFinite(u) || !isFinite(v)) {
        return Vector{0.0F};
    }
    if (!isFinite(depth)) {
        depth = 0.0F;
    }

    depth = std::clamp(depth, 0.0F, static_cast<float>(mipLevelCount - 1));
    const auto level = static_cast<std::size_t>(depth);
    const std::size_t nextLevel = std::min(level + 1, mipLevelCount - 1);
    const float levelBlend = depth - static_cast<float>(level);
    v = 1.0F - v;

    const int levelWidth = mipExtent(width, level);
    const int levelHeight = mipExtent(height, level);
    const int nextWidth = mipExtent(width, nextLevel);
    const int nextHeight = mipExtent(height, nextLevel);
    const Vector first = bilinearSample<Vector>(levels[level], levelWidth, levelHeight, channels, u, v);
    const Vector second =
        bilinearSample<Vector>(levels[nextLevel], nextWidth, nextHeight, channels, u, v);
    return first * (1.0F - levelBlend) + second * levelBlend;
}

}  // namespace

vec3 ImageData::sampleRgb(float u, float v, float depth) const {
    return sampleImage<vec3>(levels_, mipLevelCount_, width_, height_, channels_, u, v, depth);
}

vec4 ImageData::sampleRgba(float u, float v, float depth) const {
    return sampleImage<vec4>(levels_, mipLevelCount_, width_, height_, channels_, u, v, depth);
}

bool ImageData::empty() const noexcept {
    return mipLevelCount_ == 0 || levels_[0].empty();
}

bool ImageData::hasTransparency() const noexcept {
    if (empty() || (channels_ != 2 && channels_ != 4)) {
        return false;
    }
    const int alphaChannel = channels_ - 1;
    for (std::size_t index = static_cast<std::size_t>(alphaChannel); index < levels_[0].size();
         index += static_cast<std::size_t>(channels_)) {
        if (levels_[0][index] < 255U) {
            return true;
        }
    }
    return false;
}

int ImageData::width() const noexcept {
    return width_;
}

int ImageData::height() const noexcept {
    return height_;
}

int ImageData::channels() const noexcept {
    return channels_;
}

std::size_t ImageData::mipLevelCount() const noexcept {
    return mipLevelCount_;
}

int ImageData::mipWidth(std::size_t level) const {
    if (level >= mipLevelCount_) {
        throw std::out_of_range("Image mip level is out of range.");
    }
    return mipExtent(width_, level);
}

int ImageData::mipHeight(std::size_t level) const {
    if (level >= mipLevelCount_) {
        throw std::out_of_range("Image mip level is out of range.");
    }
    return mipExtent(height_, level);
}

const std::vector<std::uint8_t>& ImageData::mipPixels(std::size_t level) const {
    if (level >= mipLevelCount_) {
        throw std::out_of_range("Image mip level is out of range.");
    }
    return levels_[level];
}

std::string urlDecode(const std::string& source) {
    std::string decoded;
    decoded.reserve(source.size());
    for (std::size_t index = 0; index < source.size(); ++index) {
        if (source[index] == '%' && index + 2 < source.size() && isHexDigit(source[index + 1]) &&
            isHexDigit(source[index + 2])) {
            const unsigned int value = hexValue(source[index + 1]) * 16U + hexValue(source[index + 2]);
            decoded.push_back(static_cast<char>(value));
            index += 2;
        } else {
            decoded.push_back(source[index]);
        }
    }
    return decoded;
}

ImageData Material::loadDds(const std::filesystem::path& filename) {
    detail::ByteImage decoded = detail::loadByteImage(filename, "dds_to_array");
    ImageData image;
    image.setBaseLevel(
        decoded.width, decoded.height, decoded.channels, std::move(decoded.pixels));
    image.generateMipmaps();
    return image;
}

ImageData Material::loadPng(const std::filesystem::path& filename) {
    FileHandle file{filename};
    PngImageHandle handle;
    png_image& png = handle.get();
    if (png_image_begin_read_from_stdio(&png, file.get()) == 0) {
        const std::string error = png.message;
        throw std::runtime_error("Failed to read PNG header from " + filename.string() + ": " + error);
    }

    png.format = PNG_FORMAT_RGBA;
    const std::size_t width = static_cast<std::size_t>(png.width);
    const std::size_t height = static_cast<std::size_t>(png.height);
    constexpr std::size_t channelCount = 4;
    const std::size_t maximumSize = std::numeric_limits<std::size_t>::max();
    if (width == 0 || height == 0 ||
        width > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        height > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        width > maximumSize / height || width * height > maximumSize / channelCount) {
        throw std::runtime_error("PNG dimensions are not supported: " + filename.string());
    }
    std::vector<std::uint8_t> pixels(width * height * channelCount);
    if (png_image_finish_read(&png, nullptr, pixels.data(), 0, nullptr) == 0) {
        const std::string error = png.message;
        throw std::runtime_error("Failed to decode PNG " + filename.string() + ": " + error);
    }

    ImageData image;
    image.setBaseLevel(static_cast<int>(png.width), static_cast<int>(png.height), 4, std::move(pixels));
    image.generateMipmaps();
    return image;
}

void Material::loadTexture(TextureSlot slot, const std::filesystem::path& filename) {
    if (static_cast<std::size_t>(slot) >= textures_.size()) {
        throw std::invalid_argument("Texture slot is out of range.");
    }
    std::string extension = filename.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char value) {
        return static_cast<char>(std::tolower(value));
    });

    ImageData loaded;
    if (extension == ".png") {
        loaded = loadPng(filename);
    } else if (extension == ".dds") {
        loaded = loadDds(filename);
    } else {
        throw std::runtime_error("Unsupported texture format: " + filename.string());
    }
    const std::size_t index = static_cast<std::size_t>(slot);
    textureSourcePaths_[index] = filename.lexically_normal();
    textures_[index] = std::move(loaded);
    if (slot == TextureSlot::diffuse) {
        hasCutoutTransparency = texture(slot).hasTransparency();
    }
}

bool Material::hasTexture(TextureSlot slot) const noexcept {
    const std::size_t index = static_cast<std::size_t>(slot);
    return index < textures_.size() && !textures_[index].empty();
}

const ImageData& Material::textureData(TextureSlot slot) const {
    return texture(slot);
}

const std::filesystem::path& Material::textureSourcePath(TextureSlot slot) const {
    return textureSourcePaths_.at(static_cast<std::size_t>(slot));
}

bool Material::isEmissive() const noexcept {
    const float luminance = glm::dot(glm::max(emissiveFactor, vec3{0.0F}), kLuminanceWeights);
    return hasTexture(TextureSlot::emissive) || (isFinite(luminance) && luminance > 0.0F);
}

vec4 Material::diffuseColor(float u, float v, float footprint) const {
    const ImageData& image = texture(TextureSlot::diffuse);
    if (image.empty()) {
        return vec4{glm::max(diffuseFactor, vec3{0.0F}), 1.0F};
    }
    vec4 color = image.sampleRgba(u, v, mipDepth(footprint, image.width()));
    convertToLinear(color);
    color.r *= diffuseFactor.r;
    color.g *= diffuseFactor.g;
    color.b *= diffuseFactor.b;
    return color;
}

vec3 Material::normal(float u, float v, float footprint) const {
    const ImageData& image = texture(TextureSlot::normal);
    if (image.empty()) {
        return vec3{0.0F};
    }
    return image.sampleRgb(u, v, mipDepth(footprint, image.width())) * 2.0F - 1.0F;
}

vec3 Material::emissiveColor(float u, float v, float footprint) const {
    const ImageData& image = texture(TextureSlot::emissive);
    if (image.empty()) {
        return glm::max(emissiveFactor, vec3{0.0F});
    }
    vec4 color = image.sampleRgba(u, v, mipDepth(footprint, image.width()));
    convertToLinear(color);
    return vec3{color} * glm::max(emissiveFactor, vec3{0.0F});
}

void Material::surfaceParameters(
    float u, float v, float& surfaceRoughness, float& surfaceMetallic) const {
    const ImageData& image = texture(TextureSlot::specular);
    if (image.empty()) {
        surfaceMetallic = metallic;
        surfaceRoughness = roughness;
        return;
    }
    const vec4 surfaceData = image.sampleRgba(u, v, 0.0F);
    surfaceMetallic = std::clamp(surfaceData.b, 0.0F, 0.99F);
    surfaceRoughness = std::clamp(surfaceData.g, 1.0e-3F, 1.0F);
}

const ImageData& Material::texture(TextureSlot slot) const {
    return textures_.at(static_cast<std::size_t>(slot));
}

ImageData& Material::texture(TextureSlot slot) {
    return textures_.at(static_cast<std::size_t>(slot));
}

}  // namespace raym0nade
