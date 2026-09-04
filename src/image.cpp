#include "raym0nade/image.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <png.h>

namespace raym0nade {
namespace {

constexpr float kSpatialClampThreshold = 36.0F;
constexpr float kMinimumFilterWeight = 1.0e-12F;
constexpr float kGamma = 2.2F;

[[nodiscard]] std::size_t checkedPixelCount(const int width, const int height) {
    if (width <= 0 || height <= 0) {
        throw std::invalid_argument("Film dimensions must be positive");
    }

    const auto unsignedWidth = static_cast<std::size_t>(width);
    const auto unsignedHeight = static_cast<std::size_t>(height);
    if (unsignedWidth > std::numeric_limits<std::size_t>::max() / unsignedHeight) {
        throw std::length_error("Film dimensions exceed the addressable image size");
    }
    return unsignedWidth * unsignedHeight;
}

[[nodiscard]] std::size_t pixelIndex(const int x, const int y, const int width) noexcept {
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
           static_cast<std::size_t>(x);
}

template <typename T>
void requireSize(const std::vector<T>& values, const std::size_t expected, const char* name) {
    if (values.size() != expected) {
        throw std::logic_error(std::string{name} + " does not match the film dimensions");
    }
}

[[nodiscard]] vec3 finiteOrBlack(vec3 value) noexcept {
    for (int channel = 0; channel < 3; ++channel) {
        if (!isFinite(value[channel])) {
            value[channel] = 0.0F;
        }
    }
    return value;
}

[[nodiscard]] float luminance(const vec3& color) noexcept {
    return glm::dot(color, kLuminanceWeights);
}

void spatialClampRadiance(std::vector<RadianceData>& radiance, const int width, const int height) {
    const std::size_t count = checkedPixelCount(width, height);
    requireSize(radiance, count, "Radiance buffer");

    constexpr int filterRadius = 3;
    constexpr std::array<float, 7> kernel{
        0.03125F, 0.109375F, 0.21875F, 0.28125F, 0.21875F, 0.109375F, 0.03125F};

    std::vector<float> sourceLuminance(count, 0.0F);
    std::vector<float> horizontal(count, 0.0F);
    std::vector<float> blurred(count, 0.0F);

    for (std::size_t index = 0; index < count; ++index) {
        radiance[index].radiance = finiteOrBlack(radiance[index].radiance);
        if (!isFinite(radiance[index].varianceAccumulator) ||
            radiance[index].varianceAccumulator < 0.0F) {
            radiance[index].varianceAccumulator = 0.0F;
        }
        const float value = luminance(radiance[index].radiance);
        sourceLuminance[index] = isFinite(value) ? value : 0.0F;
    }

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float sum = 0.0F;
            for (int dx = -filterRadius; dx <= filterRadius; ++dx) {
                const int neighborX = x + dx;
                if (neighborX >= 0 && neighborX < width) {
                    sum += sourceLuminance[pixelIndex(neighborX, y, width)] *
                           kernel[static_cast<std::size_t>(dx + filterRadius)];
                }
            }
            horizontal[pixelIndex(x, y, width)] = sum;
        }
    }

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float sum = 0.0F;
            for (int dy = -filterRadius; dy <= filterRadius; ++dy) {
                const int neighborY = y + dy;
                if (neighborY >= 0 && neighborY < height) {
                    sum += horizontal[pixelIndex(x, neighborY, width)] *
                           kernel[static_cast<std::size_t>(dy + filterRadius)];
                }
            }
            blurred[pixelIndex(x, y, width)] = sum;
        }
    }

    constexpr float centerContribution = kernel[filterRadius] * kernel[filterRadius];
    constexpr float nonCenterContribution = 1.0F - centerContribution;
    for (std::size_t index = 0; index < count; ++index) {
        const float value = sourceLuminance[index];
        const float neighboringValue = blurred[index] - centerContribution * value;
        if (!isFinite(value) || !isFinite(neighboringValue) || value <= kRayEpsilon ||
            neighboringValue < 0.0F || value <= kSpatialClampThreshold * neighboringValue) {
            continue;
        }

        const float factor = neighboringValue / (value + kRayEpsilon) / nonCenterContribution;
        if (!isFinite(factor) || factor < 0.0F) {
            continue;
        }

        radiance[index].radiance *= factor;
        radiance[index].varianceAccumulator *= factor * factor;
    }
}

void filterVariance(std::vector<RadianceData>& radiance, const int width, const int height) {
    constexpr int filterRadius = 1;
    constexpr std::array<float, 3> kernel{0.25F, 0.5F, 0.25F};
    const std::size_t count = checkedPixelCount(width, height);
    requireSize(radiance, count, "Radiance buffer");

    std::vector<float> filteredVariance(count, 0.0F);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            vec3 expectedRadiance{0.0F};
            float expectedSquaredRadiance = 0.0F;
            float expectedVariance = 0.0F;

            for (int dy = -filterRadius; dy <= filterRadius; ++dy) {
                for (int dx = -filterRadius; dx <= filterRadius; ++dx) {
                    const int neighborX = x + dx;
                    const int neighborY = y + dy;
                    if (neighborX < 0 || neighborX >= width || neighborY < 0 || neighborY >= height) {
                        continue;
                    }

                    const auto& neighbor = radiance[pixelIndex(neighborX, neighborY, width)];
                    if (!isFinite(neighbor.radiance) || !isFinite(neighbor.varianceAccumulator)) {
                        continue;
                    }

                    const float weight = kernel[static_cast<std::size_t>(dx + filterRadius)] *
                                         kernel[static_cast<std::size_t>(dy + filterRadius)];
                    expectedRadiance += neighbor.radiance * weight;
                    expectedSquaredRadiance += glm::dot(neighbor.radiance, neighbor.radiance) * weight;
                    expectedVariance += std::max(neighbor.varianceAccumulator, 0.0F) * weight;
                }
            }

            const float result = expectedVariance + expectedSquaredRadiance -
                                 glm::dot(expectedRadiance, expectedRadiance);
            filteredVariance[pixelIndex(x, y, width)] =
                isFinite(result) ? std::max(result, 0.0F) : 0.0F;
        }
    }

    for (std::size_t index = 0; index < count; ++index) {
        radiance[index].varianceAccumulator = filteredVariance[index];
    }
}

[[nodiscard]] float filterWeight(
    const HitInfo& centerGeometry,
    const HitInfo& neighborGeometry,
    const RadianceData& centerRadiance,
    const RadianceData& neighborRadiance,
    const bool specular) noexcept {
    constexpr float depthSigma = 1.0F;
    constexpr float normalExponent = 1024.0F;
    constexpr float radianceSigma = 1.0F;
    constexpr float materialSigma = 1.0F;
    constexpr float radianceEpsilon = 1.0e-2F;
    constexpr float roughnessEpsilon = 4.0e-2F;

    if (!isFinite(centerGeometry.position) || !isFinite(neighborGeometry.position) ||
        !isFinite(centerGeometry.surfaceNormal) || !isFinite(neighborGeometry.surfaceNormal) ||
        !isFinite(centerGeometry.shapeNormal) || !isFinite(centerRadiance.radiance) ||
        !isFinite(neighborRadiance.radiance)) {
        return 0.0F;
    }

    const float normalAlignment = std::clamp(
        glm::dot(centerGeometry.surfaceNormal, neighborGeometry.surfaceNormal), 0.0F, 1.0F);
    float weight = std::pow(normalAlignment, normalExponent);
    if (!isFinite(weight) || weight < 1.0e-6F) {
        return 0.0F;
    }

    const vec3 positionDelta = neighborGeometry.position - centerGeometry.position;
    const vec3 direction = safeNormalize(positionDelta);
    const float projected = std::clamp(
        std::abs(glm::dot(centerGeometry.shapeNormal, direction)), 0.0F, 1.0F);
    const float tangent = projected / (safeSqrt(1.0F - projected * projected) + kRayEpsilon);

    const float centerVariance = isFinite(centerRadiance.varianceAccumulator)
                                     ? std::max(centerRadiance.varianceAccumulator, 0.0F)
                                     : 0.0F;
    const float radianceDifference = glm::length(
        centerRadiance.radiance - neighborRadiance.radiance) /
        (radianceSigma * std::sqrt(centerVariance) + radianceEpsilon);

    const vec3 materialDifference{
        centerGeometry.metallic - neighborGeometry.metallic,
        centerGeometry.specular - neighborGeometry.specular,
        centerGeometry.opacity - neighborGeometry.opacity};
    const float materialDistance = glm::length(materialDifference);
    if (!isFinite(tangent) || !isFinite(radianceDifference) || !isFinite(materialDistance)) {
        return 0.0F;
    }

    const float exponent = -tangent / depthSigma - radianceDifference - materialDistance / materialSigma;
    if (exponent < -7.5F) {
        return 0.0F;
    }

    weight *= std::exp(exponent);
    if (specular) {
        weight *= std::max(centerGeometry.roughness, roughnessEpsilon);
    }
    return isFinite(weight) && weight >= 0.0F ? weight : 0.0F;
}

void filterRadianceStep(
    std::vector<RadianceData>& radiance,
    const std::vector<HitInfo>& geometry,
    const int width,
    const int height,
    const bool specular,
    const int step) {
    constexpr int filterRadius = 2;
    constexpr std::array<float, 5> kernel{0.0625F, 0.25F, 0.375F, 0.25F, 0.0625F};
    const std::size_t count = checkedPixelCount(width, height);
    requireSize(radiance, count, "Radiance buffer");
    requireSize(geometry, count, "Geometry buffer");

    std::vector<vec3> filteredRadiance(count, vec3{0.0F});
    std::vector<float> filteredVariance(count, 0.0F);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const std::size_t centerIndex = pixelIndex(x, y, width);
            const auto& centerSample = radiance[centerIndex];
            const auto& centerGeometry = geometry[centerIndex];
            if (!isFinite(centerGeometry.position)) {
                filteredRadiance[centerIndex] = finiteOrBlack(centerSample.radiance);
                filteredVariance[centerIndex] =
                    isFinite(centerSample.varianceAccumulator)
                        ? std::max(centerSample.varianceAccumulator, 0.0F)
                        : 0.0F;
                continue;
            }

            float weightSum = 0.0F;
            for (int dy = -filterRadius; dy <= filterRadius; ++dy) {
                for (int dx = -filterRadius; dx <= filterRadius; ++dx) {
                    const int neighborX = x + dx * step;
                    const int neighborY = y + dy * step;
                    if (neighborX < 0 || neighborX >= width || neighborY < 0 || neighborY >= height) {
                        continue;
                    }

                    const std::size_t neighborIndex = pixelIndex(neighborX, neighborY, width);
                    const auto& neighborSample = radiance[neighborIndex];
                    if (!isFinite(neighborSample.radiance) ||
                        !isFinite(neighborSample.varianceAccumulator)) {
                        continue;
                    }

                    float weight = kernel[static_cast<std::size_t>(dx + filterRadius)] *
                                   kernel[static_cast<std::size_t>(dy + filterRadius)];
                    if (dx != 0 || dy != 0) {
                        weight *= filterWeight(
                            centerGeometry,
                            geometry[neighborIndex],
                            centerSample,
                            neighborSample,
                            specular);
                    }
                    if (!isFinite(weight) || weight <= 0.0F) {
                        continue;
                    }

                    weightSum += weight;
                    filteredRadiance[centerIndex] += neighborSample.radiance * weight;
                    filteredVariance[centerIndex] +=
                        std::max(neighborSample.varianceAccumulator, 0.0F) * weight * weight;
                }
            }

            if (!isFinite(weightSum) || weightSum <= kMinimumFilterWeight) {
                filteredRadiance[centerIndex] = finiteOrBlack(centerSample.radiance);
                filteredVariance[centerIndex] =
                    isFinite(centerSample.varianceAccumulator)
                        ? std::max(centerSample.varianceAccumulator, 0.0F)
                        : 0.0F;
                continue;
            }

            filteredRadiance[centerIndex] /= weightSum;
            filteredVariance[centerIndex] /= weightSum * weightSum;
            if (!isFinite(filteredRadiance[centerIndex]) || !isFinite(filteredVariance[centerIndex])) {
                filteredRadiance[centerIndex] = finiteOrBlack(centerSample.radiance);
                filteredVariance[centerIndex] =
                    isFinite(centerSample.varianceAccumulator)
                        ? std::max(centerSample.varianceAccumulator, 0.0F)
                        : 0.0F;
            }
        }
    }

    for (std::size_t index = 0; index < count; ++index) {
        radiance[index].radiance = filteredRadiance[index];
        radiance[index].varianceAccumulator = filteredVariance[index];
    }
}

void filterRadiance(
    std::vector<RadianceData>& radiance,
    const std::vector<HitInfo>& geometry,
    const int width,
    const int height,
    const bool specular) {
    filterVariance(radiance, width, height);
    for (const int step : {1, 2, 4, 8, 16}) {
        filterRadianceStep(radiance, geometry, width, height, specular, step);
    }
}

struct FileCloser {
    void operator()(std::FILE* file) const noexcept {
        if (file != nullptr) {
            std::fclose(file);
        }
    }
};

using FileHandle = std::unique_ptr<std::FILE, FileCloser>;

[[nodiscard]] FileHandle openImageFile(const std::filesystem::path& filename, const bool write) {
#if defined(_WIN32) && defined(_MSC_VER)
    std::FILE* file = nullptr;
    if (_wfopen_s(&file, filename.c_str(), write ? L"wb" : L"rb") != 0) {
        file = nullptr;
    }
#elif defined(_WIN32)
    std::FILE* file = _wfopen(filename.c_str(), write ? L"wb" : L"rb");
#else
    std::FILE* file = std::fopen(filename.c_str(), write ? "wb" : "rb");
#endif
    return FileHandle{file};
}

class PngImageHandle {
public:
    PngImageHandle() noexcept {
        image.version = PNG_IMAGE_VERSION;
    }

    PngImageHandle(const PngImageHandle&) = delete;
    PngImageHandle& operator=(const PngImageHandle&) = delete;

    ~PngImageHandle() {
        png_image_free(&image);
    }

    png_image image{};
};

[[nodiscard]] std::string pngError(const char* operation, const png_image& image) {
    const char* message = image.message[0] == '\0' ? "unknown libpng error" : image.message;
    return std::string{operation} + ": " + message;
}

void accumulateRadiance(RadianceData& destination, const vec3& incomingRadiance, const float weight) noexcept {
    if (!isFinite(incomingRadiance) || !isFinite(weight) || weight <= 0.0F) {
        return;
    }

    const vec3 weightedRadiance = incomingRadiance * weight;
    const vec3 accumulatedRadiance = destination.radiance + weightedRadiance;
    if (!isFinite(weightedRadiance) || !isFinite(accumulatedRadiance)) {
        return;
    }

    const double squaredMagnitude =
        static_cast<double>(incomingRadiance.x) * incomingRadiance.x +
        static_cast<double>(incomingRadiance.y) * incomingRadiance.y +
        static_cast<double>(incomingRadiance.z) * incomingRadiance.z;
    const double accumulatedVariance = static_cast<double>(destination.varianceAccumulator) +
                                       squaredMagnitude * static_cast<double>(weight);
    if (!std::isfinite(accumulatedVariance) ||
        accumulatedVariance > static_cast<double>(std::numeric_limits<float>::max())) {
        return;
    }

    destination.radiance = accumulatedRadiance;
    destination.varianceAccumulator = static_cast<float>(accumulatedVariance);
}

}  // namespace

Film::Film(const int width, const int height) : width_(width), height_(height) {
    const std::size_t count = checkedPixelCount(width_, height_);
    gBuffer.resize(count);
    directDiffuseRadiance.resize(count);
    directSpecularRadiance.resize(count);
    indirectDiffuseRadiance.resize(count);
    indirectSpecularRadiance.resize(count);
    pixels.resize(count, vec3{0.0F});
}

Film::Film(const std::filesystem::path& filename) {
    load(filename);
}

int Film::width() const noexcept {
    return width_;
}

int Film::height() const noexcept {
    return height_;
}

void Film::spatialClamp() {
    spatialClampRadiance(directDiffuseRadiance, width_, height_);
    spatialClampRadiance(directSpecularRadiance, width_, height_);
    spatialClampRadiance(indirectDiffuseRadiance, width_, height_);
    spatialClampRadiance(indirectSpecularRadiance, width_, height_);
}

void Film::filter() {
    const std::size_t count = checkedPixelCount(width_, height_);
    requireSize(gBuffer, count, "Geometry buffer");
    requireSize(directDiffuseRadiance, count, "Direct diffuse buffer");
    requireSize(directSpecularRadiance, count, "Direct specular buffer");
    requireSize(indirectDiffuseRadiance, count, "Indirect diffuse buffer");
    requireSize(indirectSpecularRadiance, count, "Indirect specular buffer");

    std::array<std::exception_ptr, 4> errors{};
    std::array<std::thread, 4> workers;
    auto launch = [&](const std::size_t index, std::vector<RadianceData>& radiance, const bool specular) {
        workers[index] = std::thread{[&, index, specular, buffer = &radiance] {
            try {
                filterRadiance(*buffer, gBuffer, width_, height_, specular);
            } catch (...) {
                errors[index] = std::current_exception();
            }
        }};
    };

    try {
        launch(0, directDiffuseRadiance, false);
        launch(1, directSpecularRadiance, true);
        launch(2, indirectDiffuseRadiance, false);
        launch(3, indirectSpecularRadiance, true);
    } catch (...) {
        for (auto& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
        throw;
    }

    for (auto& worker : workers) {
        worker.join();
    }
    for (const auto& error : errors) {
        if (error != nullptr) {
            std::rethrow_exception(error);
        }
    }
}

void Film::shade(const int options) {
    const std::size_t count = checkedPixelCount(width_, height_);
    requireSize(gBuffer, count, "Geometry buffer");
    requireSize(pixels, count, "Pixel buffer");
    requireSize(directDiffuseRadiance, count, "Direct diffuse buffer");
    requireSize(directSpecularRadiance, count, "Direct specular buffer");
    requireSize(indirectDiffuseRadiance, count, "Indirect diffuse buffer");
    requireSize(indirectSpecularRadiance, count, "Indirect specular buffer");

    for (std::size_t index = 0; index < count; ++index) {
        const HitInfo& geometry = gBuffer[index];
        if ((options & shapeNormal) != 0) {
            pixels[index] = isFinite(geometry.shapeNormal)
                                ? (geometry.shapeNormal + vec3{1.0F}) * 0.5F
                                : vec3{0.0F};
            continue;
        }
        if ((options & surfaceNormal) != 0) {
            pixels[index] = isFinite(geometry.surfaceNormal)
                                ? (geometry.surfaceNormal + vec3{1.0F}) * 0.5F
                                : vec3{0.0F};
            continue;
        }

        vec3 diffuseRadiance{0.0F};
        vec3 specularRadiance{0.0F};
        if ((options & directLight) != 0) {
            if ((options & diffuse) != 0) {
                diffuseRadiance += directDiffuseRadiance[index].radiance;
            }
            if ((options & specular) != 0) {
                specularRadiance += directSpecularRadiance[index].radiance;
            }
        }
        if ((options & indirectLight) != 0) {
            if ((options & diffuse) != 0) {
                diffuseRadiance += indirectDiffuseRadiance[index].radiance;
            }
            if ((options & specular) != 0) {
                specularRadiance += indirectSpecularRadiance[index].radiance;
            }
        }
        if ((options & (directLight | indirectLight)) == 0) {
            diffuseRadiance = vec3{1.0F};
        }

        const vec3 color = (options & baseColor) != 0 ? geometry.baseColor : vec3{1.0F};
        pixels[index] = color * diffuseRadiance + specularRadiance;
        if ((options & emission) != 0) {
            pixels[index] += geometry.emission * exposure;
        }
        pixels[index] = finiteOrBlack(pixels[index]);
    }
}

void Film::bloom() {
    constexpr int filterRadius = 2;
    constexpr std::array<float, 5> kernel{0.0625F, 0.25F, 0.375F, 0.25F, 0.0625F};
    constexpr float bloomExponent = 0.65F;
    const std::size_t count = checkedPixelCount(width_, height_);
    requireSize(pixels, count, "Pixel buffer");

    std::vector<vec3> bloomColor(count, vec3{0.0F});
    for (std::size_t index = 0; index < count; ++index) {
        const vec3 color = finiteOrBlack(pixels[index]);
        const float value = luminance(color);
        if (!isFinite(value) || value < 1.0F) {
            continue;
        }

        const float divisor = std::pow(value, bloomExponent);
        if (!isFinite(divisor) || divisor <= 0.0F) {
            continue;
        }
        bloomColor[index] = color / divisor;
        bloomColor[index] = glm::max(bloomColor[index] - vec3{1.0F}, vec3{0.0F});
    }

    for (int step = 1; step <= 16; step *= 2) {
        const std::vector<vec3> source = bloomColor;
        for (int y = 0; y < height_; ++y) {
            for (int x = 0; x < width_; ++x) {
                vec3 filtered{0.0F};
                for (int kernelY = -filterRadius; kernelY <= filterRadius; ++kernelY) {
                    for (int kernelX = -filterRadius; kernelX <= filterRadius; ++kernelX) {
                        const int neighborX = x + kernelX * step;
                        const int neighborY = y + kernelY * step;
                        if (neighborX < 0 || neighborX >= width_ || neighborY < 0 || neighborY >= height_) {
                            continue;
                        }
                        filtered += source[pixelIndex(neighborX, neighborY, width_)] *
                                    kernel[static_cast<std::size_t>(kernelY + filterRadius)] *
                                    kernel[static_cast<std::size_t>(kernelX + filterRadius)];
                    }
                }
                bloomColor[pixelIndex(x, y, width_)] = finiteOrBlack(filtered);
            }
        }
        for (std::size_t index = 0; index < count; ++index) {
            pixels[index] = finiteOrBlack(pixels[index] + bloomColor[index] / 6.0F);
        }
    }
}

void Film::depthOfFieldBlur() {
    struct PixelDepth {
        int x;
        int y;
        float depth;
    };

    constexpr float maximumCircleOfConfusion = 96.0F;
    constexpr float gainThreshold = 0.99F;
    const std::size_t count = checkedPixelCount(width_, height_);
    requireSize(gBuffer, count, "Geometry buffer");
    requireSize(pixels, count, "Pixel buffer");

    std::vector<PixelDepth> sortedPixels;
    std::vector<unsigned char> validDepth(count, 0);
    std::vector<vec3> blurredPixels(count, vec3{0.0F});
    std::vector<float> gained(count, 0.0F);
    sortedPixels.reserve(count);

    for (int y = 0; y < height_; ++y) {
        for (int x = 0; x < width_; ++x) {
            const std::size_t index = pixelIndex(x, y, width_);
            if (!isFinite(gBuffer[index].position)) {
                continue;
            }
            const float depth = glm::length(gBuffer[index].position - cameraPosition);
            if (!isFinite(depth) || depth <= kRayEpsilon) {
                continue;
            }
            validDepth[index] = 1;
            sortedPixels.push_back({x, y, depth});
        }
    }

    std::stable_sort(
        sortedPixels.begin(), sortedPixels.end(),
        [](const PixelDepth& lhs, const PixelDepth& rhs) { return lhs.depth < rhs.depth; });

    const float configuredCircle =
        isFinite(circleOfConfusion) ? std::max(circleOfConfusion, 0.0F) : 0.0F;
    const float configuredFocus = isFinite(focusDistance) ? focusDistance : 0.0F;

    for (const PixelDepth& pixel : sortedPixels) {
        const std::size_t sourceIndex = pixelIndex(pixel.x, pixel.y, width_);
        const float circle = std::clamp(
            configuredCircle * std::abs(1.0F - configuredFocus / pixel.depth) + kRayEpsilon,
            kRayEpsilon,
            maximumCircleOfConfusion);
        if (!isFinite(circle)) {
            continue;
        }

        const int radius = static_cast<int>(circle);
        double totalWeight = 0.0;
        for (int dy = -radius; dy <= radius; ++dy) {
            const float horizontalSquared =
                std::max(circle * circle - static_cast<float>(dy * dy), 0.0F);
            const int horizontalRadius = static_cast<int>(std::sqrt(horizontalSquared));

            int left = -horizontalRadius;
            for (; left <= 0; ++left) {
                const float distance = std::sqrt(static_cast<float>(left * left + dy * dy));
                const float weight = std::clamp(circle - distance, 0.0F, 1.0F);
                if (weight >= 1.0F) {
                    break;
                }
                totalWeight += weight;
            }

            int right = horizontalRadius;
            for (; right > 0; --right) {
                const float distance = std::sqrt(static_cast<float>(right * right + dy * dy));
                const float weight = std::clamp(circle - distance, 0.0F, 1.0F);
                if (weight >= 1.0F) {
                    break;
                }
                totalWeight += weight;
            }
            if (right >= left) {
                totalWeight += static_cast<double>(right - left + 1);
            }
        }

        if (!std::isfinite(totalWeight) || totalWeight <= kMinimumFilterWeight) {
            blurredPixels[sourceIndex] += finiteOrBlack(pixels[sourceIndex]);
            gained[sourceIndex] = 1.0F;
            continue;
        }

        const vec3 color = finiteOrBlack(pixels[sourceIndex]);
        for (int dy = -radius; dy <= radius; ++dy) {
            const float horizontalSquared =
                std::max(circle * circle - static_cast<float>(dy * dy), 0.0F);
            const int horizontalRadius = static_cast<int>(std::sqrt(horizontalSquared));
            for (int dx = -horizontalRadius; dx <= horizontalRadius; ++dx) {
                const int targetX = pixel.x + dx;
                const int targetY = pixel.y + dy;
                if (targetX < 0 || targetX >= width_ || targetY < 0 || targetY >= height_) {
                    continue;
                }

                const std::size_t targetIndex = pixelIndex(targetX, targetY, width_);
                if (validDepth[targetIndex] == 0) {
                    continue;
                }
                const float distance = std::sqrt(static_cast<float>(dx * dx + dy * dy));
                float weight = static_cast<float>(
                    static_cast<double>(std::clamp(circle - distance, 0.0F, 1.0F)) / totalWeight);
                if (!isFinite(weight) || weight < kRayEpsilon) {
                    continue;
                }
                if (gained[targetIndex] + weight > gainThreshold) {
                    weight = gainThreshold - gained[targetIndex];
                }
                if (weight < kRayEpsilon) {
                    continue;
                }

                blurredPixels[targetIndex] += color * weight;
                gained[targetIndex] += weight;
            }
        }
    }

    for (std::size_t index = 0; index < count; ++index) {
        if (validDepth[index] == 0 || !isFinite(gained[index]) || gained[index] <= kRayEpsilon) {
            pixels[index] = finiteOrBlack(pixels[index]);
            continue;
        }
        const vec3 result = blurredPixels[index] / gained[index];
        pixels[index] = isFinite(result) ? result : finiteOrBlack(pixels[index]);
    }
}

void Film::applyFxaa() {
    constexpr float edgeThresholdMinimum = 0.0312F;
    constexpr float edgeThresholdMaximum = 0.125F;
    constexpr float subpixelQuality = 0.75F;
    constexpr int quality = 12;
    const std::size_t count = checkedPixelCount(width_, height_);
    requireSize(pixels, count, "Pixel buffer");

    const std::vector<vec3> source = pixels;
    std::vector<vec3> output(count, vec3{0.0F});
    const auto sampleLuminance = [&](const int x, const int y) noexcept {
        const float value = luminance(finiteOrBlack(source[pixelIndex(x, y, width_)]));
        return isFinite(value) ? value : 0.0F;
    };

    for (int y = 0; y < height_; ++y) {
        for (int x = 0; x < width_; ++x) {
            const std::size_t index = pixelIndex(x, y, width_);
            const float center = sampleLuminance(x, y);
            const float north = y > 0 ? sampleLuminance(x, y - 1) : center;
            const float south = y + 1 < height_ ? sampleLuminance(x, y + 1) : center;
            const float east = x + 1 < width_ ? sampleLuminance(x + 1, y) : center;
            const float west = x > 0 ? sampleLuminance(x - 1, y) : center;

            const float rangeMinimum = std::min({north, south, east, west});
            const float rangeMaximum = std::max({north, south, east, west});
            const float range = rangeMaximum - rangeMinimum;
            if (!isFinite(range) ||
                range < std::max(edgeThresholdMinimum, rangeMaximum * edgeThresholdMaximum)) {
                output[index] = finiteOrBlack(source[index]);
                continue;
            }

            const float northwest = x > 0 && y > 0 ? sampleLuminance(x - 1, y - 1) : center;
            const float northeast = x + 1 < width_ && y > 0 ? sampleLuminance(x + 1, y - 1) : center;
            const float southwest = x > 0 && y + 1 < height_ ? sampleLuminance(x - 1, y + 1) : center;
            const float southeast = x + 1 < width_ && y + 1 < height_
                                        ? sampleLuminance(x + 1, y + 1)
                                        : center;

            const float horizontalEdge =
                std::abs((northwest + west + southwest) - (northeast + east + southeast)) / 3.0F;
            const float verticalEdge =
                std::abs((northwest + north + northeast) - (southwest + south + southeast)) / 3.0F;
            const bool horizontal = horizontalEdge >= verticalEdge;
            const float stepLength = horizontal ? 1.0F / static_cast<float>(width_)
                                                : 1.0F / static_cast<float>(height_);
            const float gradientStep = std::clamp(
                (horizontal ? horizontalEdge : verticalEdge) / range, -2.0F, 2.0F);

            const vec2 uv{
                static_cast<float>(x) / static_cast<float>(width_),
                static_cast<float>(y) / static_cast<float>(height_)};
            vec3 finalColor = finiteOrBlack(source[index]);
            float bestDifference = 0.0F;

            for (int sampleIndex = 0; sampleIndex < quality; ++sampleIndex) {
                const float distance = gradientStep * stepLength * static_cast<float>(sampleIndex + 1);
                const vec2 offset = horizontal ? vec2{0.0F, distance} : vec2{distance, 0.0F};
                const vec2 sampleUv = uv + offset;
                if (sampleUv.x < 0.0F || sampleUv.x > 1.0F ||
                    sampleUv.y < 0.0F || sampleUv.y > 1.0F) {
                    continue;
                }

                const int sampleX = std::clamp(
                    static_cast<int>(sampleUv.x * static_cast<float>(width_)), 0, width_ - 1);
                const int sampleY = std::clamp(
                    static_cast<int>(sampleUv.y * static_cast<float>(height_)), 0, height_ - 1);
                const float difference = std::abs(sampleLuminance(sampleX, sampleY) - center);
                if (difference > bestDifference) {
                    bestDifference = difference;
                    finalColor = finiteOrBlack(source[pixelIndex(sampleX, sampleY, width_)]);
                }
            }

            const float subpixelOffset = std::min(
                (std::abs(north + south - 2.0F * center) * 2.0F +
                 std::abs(east + west - 2.0F * center)) *
                    0.25F,
                1.0F);
            output[index] = finiteOrBlack(glm::mix(
                finiteOrBlack(source[index]), finalColor, subpixelOffset * subpixelQuality));
        }
    }

    pixels = std::move(output);
}

void Film::gammaCorrect() {
    const std::size_t count = checkedPixelCount(width_, height_);
    requireSize(pixels, count, "Pixel buffer");

    for (vec3& pixel : pixels) {
        pixel = glm::max(finiteOrBlack(pixel), vec3{0.0F});
        const double value = 0.3 * static_cast<double>(pixel.x) +
                             0.6 * static_cast<double>(pixel.y) +
                             0.1 * static_cast<double>(pixel.z);
        if (std::isfinite(value) && value > 0.75) {
            const double boundedValue = std::tanh(3.0 * (value - 0.75)) / 3.0 + 0.75;
            pixel *= static_cast<float>(boundedValue / value);
        } else if (!std::isfinite(value)) {
            pixel = vec3{1.0F};
        }
        pixel = glm::clamp(pixel, vec3{0.0F}, vec3{1.0F});
        pixel = glm::pow(pixel, vec3{1.0F / kGamma});
    }
}

void Film::reverseGammaCorrect() {
    const std::size_t count = checkedPixelCount(width_, height_);
    requireSize(pixels, count, "Pixel buffer");
    for (vec3& pixel : pixels) {
        pixel = glm::pow(glm::clamp(finiteOrBlack(pixel), vec3{0.0F}, vec3{1.0F}), vec3{kGamma});
    }
}

void Film::postProcess(const int shadeOptions) {
    shade(shadeOptions);
    if ((shadeOptions & depthOfFieldEnabled) != 0) {
        depthOfFieldBlur();
    }
    if ((shadeOptions & bloomEnabled) != 0) {
        bloom();
    }
    gammaCorrect();
    if ((shadeOptions & fxaaEnabled) != 0) {
        applyFxaa();
    }
}

void Film::save(const std::filesystem::path& filename) const {
    const std::size_t count = checkedPixelCount(width_, height_);
    requireSize(pixels, count, "Pixel buffer");

    if (count > std::numeric_limits<std::size_t>::max() / 3U) {
        throw std::length_error("PNG output buffer is too large");
    }
    std::vector<png_byte> imageData(count * 3U, 0);
    for (std::size_t index = 0; index < count; ++index) {
        const vec3 color = glm::clamp(finiteOrBlack(pixels[index]), vec3{0.0F}, vec3{1.0F});
        for (int channel = 0; channel < 3; ++channel) {
            imageData[index * 3U + static_cast<std::size_t>(channel)] =
                static_cast<png_byte>(color[channel] * 255.0F);
        }
    }

    FileHandle file = openImageFile(filename, true);
    if (file == nullptr) {
        throw std::runtime_error("Could not open PNG file for writing: " + filename.string());
    }

    PngImageHandle png;
    png.image.width = static_cast<png_uint_32>(width_);
    png.image.height = static_cast<png_uint_32>(height_);
    png.image.format = PNG_FORMAT_RGB;
    if (png_image_write_to_stdio(&png.image, file.get(), 0, imageData.data(), 0, nullptr) == 0) {
        throw std::runtime_error(pngError("Could not write PNG image", png.image));
    }
}

void Film::load(const std::filesystem::path& filename) {
    FileHandle file = openImageFile(filename, false);
    if (file == nullptr) {
        throw std::runtime_error("Could not open PNG file for reading: " + filename.string());
    }

    PngImageHandle png;
    if (png_image_begin_read_from_stdio(&png.image, file.get()) == 0) {
        throw std::runtime_error(pngError("Could not read PNG header", png.image));
    }
    if (png.image.width > static_cast<png_uint_32>(std::numeric_limits<int>::max()) ||
        png.image.height > static_cast<png_uint_32>(std::numeric_limits<int>::max())) {
        throw std::length_error("PNG dimensions exceed the supported integer range");
    }

    const int loadedWidth = static_cast<int>(png.image.width);
    const int loadedHeight = static_cast<int>(png.image.height);
    const std::size_t count = checkedPixelCount(loadedWidth, loadedHeight);
    if (count > std::numeric_limits<std::size_t>::max() / 4U) {
        throw std::length_error("PNG input buffer is too large");
    }

    png.image.format = PNG_FORMAT_RGBA;
    std::vector<png_byte> imageData(count * 4U, 0);
    if (png_image_finish_read(&png.image, nullptr, imageData.data(), 0, nullptr) == 0) {
        throw std::runtime_error(pngError("Could not decode PNG image", png.image));
    }

    std::vector<vec3> loadedPixels(count, vec3{0.0F});
    for (std::size_t index = 0; index < count; ++index) {
        loadedPixels[index] = vec3{
            static_cast<float>(imageData[index * 4U]) / 255.0F,
            static_cast<float>(imageData[index * 4U + 1U]) / 255.0F,
            static_cast<float>(imageData[index * 4U + 2U]) / 255.0F};
    }

    width_ = loadedWidth;
    height_ = loadedHeight;
    pixels = std::move(loadedPixels);
    gBuffer.clear();
    directDiffuseRadiance.clear();
    directSpecularRadiance.clear();
    indirectDiffuseRadiance.clear();
    indirectSpecularRadiance.clear();
}

void accumulateInwardRadiance(
    const vec3& baseColor,
    const LightSample& sample,
    RadianceData& diffuseRadiance,
    RadianceData& specularRadiance) {
    if (!isFinite(baseColor) || !isFinite(sample.radiance) || !isFinite(sample.throughput) ||
        !isFinite(sample.weight) || sample.weight <= 0.0F ||
        glm::length(sample.radiance) < kRayEpsilon) {
        return;
    }

    const float baseLength = glm::length(baseColor);
    if (!isFinite(baseLength) || baseLength < kRayEpsilon) {
        accumulateRadiance(
            specularRadiance, sample.radiance * sample.throughput, sample.weight);
        return;
    }

    constexpr float inverseSqrtThree = 0.5773502691896258F;
    constexpr vec3 white{inverseSqrtThree};
    const vec3 normalizedBase = safeNormalize(baseColor);
    const float alignment = std::clamp(glm::dot(normalizedBase, white), -1.0F, 1.0F);
    const float minimumComponent = std::min({
        std::abs(baseColor.x), std::abs(baseColor.y), std::abs(baseColor.z)});
    if (alignment > 0.99F && minimumComponent > kRayEpsilon) {
        accumulateRadiance(
            diffuseRadiance,
            sample.radiance * (sample.throughput / baseColor),
            sample.weight);
        return;
    }

    const vec3 perpendicularVector = glm::cross(normalizedBase, white);
    if (glm::dot(perpendicularVector, perpendicularVector) <= kRayEpsilon * kRayEpsilon) {
        accumulateRadiance(
            specularRadiance, sample.radiance * sample.throughput, sample.weight);
        return;
    }
    const vec3 perpendicular = safeNormalize(perpendicularVector);
    const vec3 planarThroughput =
        sample.throughput - perpendicular * glm::dot(perpendicular, sample.throughput);

    const float sumDenominator = 1.0F + alignment;
    const float differenceDenominator = 1.0F - alignment;
    if (std::abs(sumDenominator) <= kRayEpsilon ||
        std::abs(differenceDenominator) <= kRayEpsilon) {
        accumulateRadiance(
            specularRadiance, sample.radiance * sample.throughput, sample.weight);
        return;
    }

    const float projectionOnWhite = glm::dot(planarThroughput, white);
    const float projectionOnBase = glm::dot(planarThroughput, normalizedBase);
    const float sum = (projectionOnWhite + projectionOnBase) / sumDenominator;
    const float difference = (projectionOnWhite - projectionOnBase) / differenceDenominator;
    const float baseCoefficient = 0.5F * (sum - difference);
    if (!isFinite(baseCoefficient)) {
        return;
    }

    const vec3 diffuseThroughput = baseCoefficient * normalizedBase;
    accumulateRadiance(
        diffuseRadiance,
        sample.radiance * (baseCoefficient / baseLength),
        sample.weight);
    accumulateRadiance(
        specularRadiance,
        sample.radiance * (sample.throughput - diffuseThroughput),
        sample.weight);
}

}  // namespace raym0nade
