#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "raym0nade/geometry.hpp"

namespace raym0nade {

struct ImageExtent {
    std::uint32_t width{640U};
    std::uint32_t height{360U};

    void validate() const;
    [[nodiscard]] std::size_t pixelCount() const;
};

struct PinholeCamera {
    vec3 position{0.0F};
    vec3 direction{0.0F, 0.0F, 1.0F};
    vec3 up{0.0F, 1.0F, 0.0F};
    vec3 right{1.0F, 0.0F, 0.0F};
    float pixelScale{1.0e-3F};

    void validate() const;
};

enum class PrimaryAov : std::uint32_t {
    BaseColor,
    ShapeNormal,
};

struct PrimaryRenderRequest {
    ImageExtent extent{};
    PinholeCamera camera{};
    PrimaryAov aov{PrimaryAov::BaseColor};

    void validate() const;
};

struct LinearImage {
    ImageExtent extent{};
    std::vector<vec3> pixels;

    void validate() const;
};

// The camera and extent must have been validated before entering a pixel loop.
// Pixel coordinates use the legacy integer sample location without a half-pixel offset.
[[nodiscard]] vec3 primaryRayDirectionUnnormalized(
    const PinholeCamera& camera,
    const ImageExtent& extent,
    std::uint32_t x,
    std::uint32_t y) noexcept;
[[nodiscard]] Ray makePrimaryRay(
    const PinholeCamera& camera,
    const ImageExtent& extent,
    std::uint32_t x,
    std::uint32_t y) noexcept;

}  // namespace raym0nade
