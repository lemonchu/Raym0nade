#include <array>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "raym0nade/render.hpp"

#ifndef RAYM0NADE_TEST_SOURCE_DIR
#error "RAYM0NADE_TEST_SOURCE_DIR must be defined by CMake."
#endif

namespace {

using namespace raym0nade;

void expect(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void expectNear(float actual, float expected, float tolerance, const std::string& message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message);
    }
}

void expectNear(
    const vec3& actual,
    const vec3& expected,
    float tolerance,
    const std::string& message) {
    expectNear(actual.x, expected.x, tolerance, message + " (x)");
    expectNear(actual.y, expected.y, tolerance, message + " (y)");
    expectNear(actual.z, expected.z, tolerance, message + " (z)");
}

template <typename Function>
void expectInvalidArgument(Function&& function, const std::string& message) {
    try {
        function();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

bool exactlyEqual(const vec3& lhs, const vec3& rhs) noexcept {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

void testContractValidation() {
    ImageExtent extent{};
    extent.width = 0U;
    expectInvalidArgument([&] { extent.validate(); }, "A zero-width image must be rejected.");

    PinholeCamera camera{};
    camera.right = camera.direction;
    expectInvalidArgument(
        [&] { camera.validate(); }, "A linearly dependent camera basis must be rejected.");

    PrimaryRenderRequest request{};
    request.aov = static_cast<PrimaryAov>(std::numeric_limits<std::uint32_t>::max());
    expectInvalidArgument(
        [&] { request.validate(); }, "An unknown primary AOV must be rejected.");

    LinearImage image{ImageExtent{2U, 2U}, std::vector<vec3>(3U, vec3{0.0F})};
    expectInvalidArgument(
        [&] { image.validate(); }, "A mismatched linear image pixel count must be rejected.");
}

void testLegacyRayConvention() {
    const ImageExtent extent{4U, 4U};
    PinholeCamera camera{};
    camera.pixelScale = 0.20F;
    camera.right = vec3{2.0F, 0.0F, 0.0F};
    camera.up = vec3{0.0F, 3.0F, 0.0F};
    camera.validate();

    const Ray center = makePrimaryRay(camera, extent, 2U, 2U);
    expectNear(center.origin, vec3{0.0F}, 0.0F, "The primary ray must preserve camera position.");
    expectNear(
        center.direction,
        vec3{0.0F, 0.0F, 1.0F},
        0.0F,
        "The integer center pixel must not receive a half-pixel offset.");

    const Ray right = makePrimaryRay(camera, extent, 3U, 2U);
    expectNear(
        right.direction,
        safeNormalize(vec3{0.4F, 0.0F, 1.0F}),
        1.0e-7F,
        "The right camera vector must retain its original magnitude.");

    const Ray lower = makePrimaryRay(camera, extent, 2U, 3U);
    expectNear(
        lower.direction,
        safeNormalize(vec3{0.0F, 0.6F, 1.0F}),
        1.0e-7F,
        "Increasing the raster row must retain the legacy positive-up convention.");
}

void testPrimaryAovs() {
    const std::filesystem::path sourceDirectory{RAYM0NADE_TEST_SOURCE_DIR};
    Model model{sourceDirectory / "tests" / "data", "triangle.obj", "null"};

    PrimaryRenderRequest request{};
    request.extent = ImageExtent{4U, 4U};
    request.camera.pixelScale = 0.20F;
    request.aov = PrimaryAov::BaseColor;

    const LinearImage baseColor = renderPrimaryAovCpu(model, request);
    baseColor.validate();
    expect(
        baseColor.extent.width == 4U && baseColor.extent.height == 4U,
        "The CPU primary renderer must preserve the requested dimensions.");
    expect(baseColor.pixels.size() == 16U, "A 4x4 primary render must contain 16 pixels.");

    constexpr std::array<const char*, 4> expectedMask{
        "1111",
        "0111",
        "0111",
        "0010",
    };
    for (std::uint32_t y = 0U; y < request.extent.height; ++y) {
        for (std::uint32_t x = 0U; x < request.extent.width; ++x) {
            const std::size_t index =
                static_cast<std::size_t>(y) * request.extent.width + x;
            const bool expectedHit = expectedMask[y][x] == '1';
            expectNear(
                baseColor.pixels[index],
                expectedHit ? vec3{0.8F, 0.2F, 0.1F} : vec3{0.0F},
                1.0e-6F,
                "The base-color AOV does not match the expected triangle mask.");
        }
    }

    const LinearImage repeated = renderPrimaryAovCpu(model, request);
    expect(repeated.pixels.size() == baseColor.pixels.size(), "Repeated render size changed.");
    for (std::size_t index = 0; index < baseColor.pixels.size(); ++index) {
        expect(
            exactlyEqual(repeated.pixels[index], baseColor.pixels[index]),
            "Repeated CPU primary renders must be pixel-identical.");
    }

    request.aov = PrimaryAov::ShapeNormal;
    const LinearImage shapeNormal = renderPrimaryAovCpu(model, request);
    shapeNormal.validate();
    for (std::uint32_t y = 0U; y < request.extent.height; ++y) {
        for (std::uint32_t x = 0U; x < request.extent.width; ++x) {
            const std::size_t index =
                static_cast<std::size_t>(y) * request.extent.width + x;
            const bool expectedHit = expectedMask[y][x] == '1';
            expectNear(
                shapeNormal.pixels[index],
                expectedHit ? vec3{0.5F, 0.5F, 0.0F} : vec3{0.0F},
                1.0e-6F,
                "The encoded shape-normal AOV does not match the expected triangle mask.");
        }
    }
}

}  // namespace

int main() {
    try {
        testContractValidation();
        testLegacyRayConvention();
        testPrimaryAovs();
        std::cout << "Primary render contract test passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
