#include "raym0nade/component.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

#include "python_image_loader.hpp"

namespace raym0nade {
namespace {

constexpr std::uint32_t kFloatRandomBits = 24U;
constexpr std::uint32_t kFloatRandomShift = std::mt19937::word_size - kFloatRandomBits;
constexpr double kOpenUnitDenominator = static_cast<double>((std::uint32_t{1} << kFloatRandomBits) + 1U);

float sanitizedRadiance(float value) noexcept {
    return std::isfinite(value) && value > 0.0F ? value : 0.0F;
}

}  // namespace

Generator::Generator(std::uint32_t seed) : engine_(seed) {}

float Generator::operator()() noexcept {
    const std::uint32_t uniformlyDistributedBits = engine_() >> kFloatRandomShift;
    return static_cast<float>(
        (static_cast<double>(uniformlyDistributedBits) + 1.0) / kOpenUnitDenominator);
}

void RandomDistribution::initialize(const std::vector<float>& weights) {
    cumulativeWeights_.clear();
    cumulativeWeights_.reserve(weights.size());
    totalWeight_ = 0.0;

    for (const float weight : weights) {
        if (std::isfinite(weight) && weight > 0.0F) {
            totalWeight_ += static_cast<double>(weight);
        }
        cumulativeWeights_.push_back(totalWeight_);
    }

    if (!std::isfinite(totalWeight_) || !(totalWeight_ > 0.0)) {
        totalWeight_ = 0.0;
        std::fill(cumulativeWeights_.begin(), cumulativeWeights_.end(), 0.0);
    }
}

int RandomDistribution::sample(Generator& generator) const noexcept {
    if (empty()) {
        return -1;
    }

    const double target = totalWeight_ * static_cast<double>(generator());
    const auto selected = std::upper_bound(cumulativeWeights_.begin(), cumulativeWeights_.end(), target);
    if (selected == cumulativeWeights_.end()) {
        return -1;
    }
    return static_cast<int>(std::distance(cumulativeWeights_.begin(), selected));
}

int RandomDistribution::operator()(Generator& generator) const noexcept {
    return sample(generator);
}

float RandomDistribution::pdf(std::size_t index) const noexcept {
    if (empty() || index >= cumulativeWeights_.size()) {
        return 0.0F;
    }

    const double previous = index == 0 ? 0.0 : cumulativeWeights_[index - 1];
    const double weight = cumulativeWeights_[index] - previous;
    return weight > 0.0 ? static_cast<float>(weight / totalWeight_) : 0.0F;
}

bool RandomDistribution::empty() const noexcept {
    return cumulativeWeights_.empty() || !(totalWeight_ > 0.0) || !std::isfinite(totalWeight_);
}

std::size_t RandomDistribution::size() const noexcept {
    return cumulativeWeights_.size();
}

vec3 Face::center() const noexcept {
    return (vertices[0] + vertices[1] + vertices[2]) / 3.0F;
}

Box Face::bounds() const noexcept {
    return Box{
        glm::min(vertices[0], glm::min(vertices[1], vertices[2])),
        glm::max(vertices[0], glm::max(vertices[1], vertices[2])),
    };
}

HitRecord::HitRecord(float minimumDistance, float maximumDistance)
    : tMinimum(minimumDistance), tMaximum(maximumDistance) {}

void SkyBox::initializeDistribution() {
    const std::size_t width = static_cast<std::size_t>(width_);
    const std::size_t height = static_cast<std::size_t>(height_);
    const std::size_t pixelCount = width * height;

    solidAngles_.assign(pixelCount, 0.0F);
    std::vector<float> weights(pixelCount, 0.0F);

    const double azimuthWidth = 2.0 * static_cast<double>(kPi) / static_cast<double>(width_);
    for (int row = 0; row < height_; ++row) {
        const double polarMinimum = static_cast<double>(kPi) * static_cast<double>(row) /
                                    static_cast<double>(height_);
        const double polarMaximum = static_cast<double>(kPi) * static_cast<double>(row + 1) /
                                    static_cast<double>(height_);
        const float solidAngle = static_cast<float>(
            azimuthWidth * (std::cos(polarMinimum) - std::cos(polarMaximum)));

        for (int column = 0; column < width_; ++column) {
            const std::size_t index = static_cast<std::size_t>(row) * width +
                                      static_cast<std::size_t>(column);
            solidAngles_[index] = solidAngle;
            const float luminance = glm::dot(radiance_[index], kLuminanceWeights);
            weights[index] = std::isfinite(luminance) && luminance > 0.0F
                                 ? luminance * solidAngle
                                 : 0.0F;
        }
    }

    distribution_.initialize(weights);
}

void SkyBox::load(const std::filesystem::path& filename) {
    const detail::FloatImage image = detail::loadFloatImage(filename, "hdr_to_array");
    if (image.width <= 0 || image.height <= 0 || image.channels < 3) {
        throw std::runtime_error("Decoded HDR image has invalid dimensions or channel count.");
    }

    const std::size_t width = static_cast<std::size_t>(image.width);
    const std::size_t height = static_cast<std::size_t>(image.height);
    if (width > std::numeric_limits<std::size_t>::max() / height) {
        throw std::runtime_error("Decoded HDR image dimensions are too large.");
    }
    const std::size_t pixelCount = width * height;
    if (pixelCount > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Decoded HDR image contains too many pixels.");
    }
    const std::size_t channelCount = static_cast<std::size_t>(image.channels);
    if (pixelCount > std::numeric_limits<std::size_t>::max() / channelCount ||
        image.pixels.size() != pixelCount * channelCount) {
        throw std::runtime_error("Decoded HDR pixel data does not match its dimensions.");
    }

    SkyBox loaded;
    loaded.width_ = image.width;
    loaded.height_ = image.height;
    loaded.radiance_.resize(pixelCount);
    for (std::size_t index = 0; index < pixelCount; ++index) {
        const std::size_t offset = index * channelCount;
        loaded.radiance_[index] = vec3{
            sanitizedRadiance(image.pixels[offset]),
            sanitizedRadiance(image.pixels[offset + 1]),
            sanitizedRadiance(image.pixels[offset + 2]),
        };
    }
    loaded.initializeDistribution();
    *this = std::move(loaded);
}

bool SkyBox::empty() const noexcept {
    return width_ <= 0 || height_ <= 0 || radiance_.empty();
}

int SkyBox::width() const noexcept {
    return width_;
}

int SkyBox::height() const noexcept {
    return height_;
}

const std::vector<vec3>& SkyBox::radiancePixels() const noexcept {
    return radiance_;
}

vec3 SkyBox::radiance(const vec3& direction) const noexcept {
    if (empty() || !isFinite(direction)) {
        return vec3{0.0F};
    }

    const vec3 unitDirection = safeNormalize(direction);
    if (glm::dot(unitDirection, unitDirection) == 0.0F) {
        return vec3{0.0F};
    }

    double azimuth = std::atan2(
        -static_cast<double>(unitDirection.x), static_cast<double>(unitDirection.z));
    if (azimuth < 0.0) {
        azimuth += 2.0 * static_cast<double>(kPi);
    }
    const double polar = std::acos(std::clamp(static_cast<double>(unitDirection.y), -1.0, 1.0));

    const int column = std::clamp(
        static_cast<int>(azimuth / (2.0 * static_cast<double>(kPi)) * static_cast<double>(width_)),
        0,
        width_ - 1);
    const int row = std::clamp(
        static_cast<int>(polar / static_cast<double>(kPi) * static_cast<double>(height_)),
        0,
        height_ - 1);
    const std::size_t index = static_cast<std::size_t>(row) * static_cast<std::size_t>(width_) +
                              static_cast<std::size_t>(column);
    return radiance_[index];
}

bool SkyBox::sample(
    Generator& generator, int& pixelIndex, vec3& direction, vec3& weightedRadiance) const noexcept {
    pixelIndex = -1;
    direction = vec3{0.0F};
    weightedRadiance = vec3{0.0F};

    const int selected = distribution_.sample(generator);
    if (selected < 0 || static_cast<std::size_t>(selected) >= radiance_.size()) {
        return false;
    }

    const float probability = distribution_.pdf(static_cast<std::size_t>(selected));
    const float solidAngle = solidAngles_[static_cast<std::size_t>(selected)];
    if (!(probability > 0.0F) || !(solidAngle > 0.0F) || !std::isfinite(probability) ||
        !std::isfinite(solidAngle)) {
        return false;
    }

    const int column = selected % width_;
    const int row = selected / width_;
    const float polar = kPi * (static_cast<float>(row) + 0.5F) / static_cast<float>(height_);
    const float azimuth = 2.0F * kPi * (static_cast<float>(column) + 0.5F) /
                          static_cast<float>(width_);

    direction = vec3{
        -std::sin(polar) * std::sin(azimuth),
        std::cos(polar),
        std::sin(polar) * std::cos(azimuth),
    };
    weightedRadiance = radiance_[static_cast<std::size_t>(selected)] * (solidAngle / probability);
    if (!isFinite(direction) || !isFinite(weightedRadiance)) {
        direction = vec3{0.0F};
        weightedRadiance = vec3{0.0F};
        return false;
    }

    pixelIndex = selected;
    return true;
}

}  // namespace raym0nade
