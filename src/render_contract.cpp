#include "raym0nade/render_contract.hpp"

#include <cmath>
#include <stdexcept>

namespace raym0nade {
namespace {

constexpr std::uint64_t kMaximumPixelCount = 1ULL << 30U;

}  // namespace

void ImageExtent::validate() const {
    if (width == 0U || height == 0U) {
        throw std::invalid_argument("Image width and height must be positive.");
    }
    const std::uint64_t count =
        static_cast<std::uint64_t>(width) * static_cast<std::uint64_t>(height);
    if (count > kMaximumPixelCount) {
        throw std::invalid_argument("Image dimensions are too large.");
    }
}

std::size_t ImageExtent::pixelCount() const {
    validate();
    return static_cast<std::size_t>(width) * static_cast<std::size_t>(height);
}

void PinholeCamera::validate() const {
    if (!isFinite(position) || !isFinite(direction) || !isFinite(up) || !isFinite(right) ||
        glm::dot(direction, direction) <= 1.0e-12F || glm::dot(up, up) <= 1.0e-12F ||
        glm::dot(right, right) <= 1.0e-12F) {
        throw std::invalid_argument("Camera vectors must be finite and non-zero.");
    }

    const vec3 normalizedDirection = safeNormalize(direction);
    const vec3 normalizedRight = safeNormalize(right);
    const vec3 normalizedUp = safeNormalize(up);
    const float basisVolume = std::abs(
        glm::dot(glm::cross(normalizedDirection, normalizedRight), normalizedUp));
    if (!isFinite(basisVolume) || basisVolume <= 1.0e-4F) {
        throw std::invalid_argument(
            "Camera direction, right, and up vectors must be linearly independent.");
    }
    if (!isFinite(pixelScale) || pixelScale <= 0.0F) {
        throw std::invalid_argument("Camera pixel scale must be finite and positive.");
    }
}

void PrimaryRenderRequest::validate() const {
    extent.validate();
    camera.validate();
    switch (aov) {
        case PrimaryAov::BaseColor:
        case PrimaryAov::ShapeNormal:
            return;
    }
    throw std::invalid_argument("Primary AOV is not supported.");
}

void LinearImage::validate() const {
    const std::size_t expectedPixelCount = extent.pixelCount();
    if (pixels.size() != expectedPixelCount) {
        throw std::invalid_argument("Linear image pixel count does not match its extent.");
    }
    for (const vec3& pixel : pixels) {
        if (!isFinite(pixel)) {
            throw std::invalid_argument("Linear image pixels must be finite.");
        }
    }
}

vec3 primaryRayDirectionUnnormalized(
    const PinholeCamera& camera,
    const ImageExtent& extent,
    std::uint32_t x,
    std::uint32_t y) noexcept {
    const float imageX = static_cast<float>(x) - static_cast<float>(extent.width) * 0.5F;
    const float imageY = static_cast<float>(y) - static_cast<float>(extent.height) * 0.5F;
    return camera.direction +
           camera.pixelScale * (imageX * camera.right + imageY * camera.up);
}

Ray makePrimaryRay(
    const PinholeCamera& camera,
    const ImageExtent& extent,
    std::uint32_t x,
    std::uint32_t y) noexcept {
    return Ray{
        camera.position,
        safeNormalize(primaryRayDirectionUnnormalized(camera, extent, x, y)),
    };
}

}  // namespace raym0nade
