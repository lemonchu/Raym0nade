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

const vec3& pixelAt(const LinearImage& image, std::uint32_t x, std::uint32_t y) {
    return image.pixels[static_cast<std::size_t>(y) * image.extent.width + x];
}

void expectImagesExactlyEqual(
    const LinearImage& actual, const LinearImage& expected, const std::string& message) {
    expect(
        actual.extent.width == expected.extent.width &&
            actual.extent.height == expected.extent.height &&
            actual.pixels.size() == expected.pixels.size(),
        message + " (dimensions)");
    for (std::size_t index = 0U; index < actual.pixels.size(); ++index) {
        expect(exactlyEqual(actual.pixels[index], expected.pixels[index]), message);
    }
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
    request.directionalLight.directionToLight = vec3{0.0F};
    expectInvalidArgument(
        [&] { request.validate(); },
        "A zero directional-light direction must be rejected for every primary AOV.");

    request = PrimaryRenderRequest{};
    request.directionalLight.directionToLight.x =
        std::numeric_limits<float>::quiet_NaN();
    expectInvalidArgument(
        [&] { request.validate(); }, "A non-finite directional-light direction must be rejected.");

    request = PrimaryRenderRequest{};
    request.directionalLight.incidentRadiance.y = -1.0F;
    expectInvalidArgument(
        [&] { request.validate(); }, "Negative directional-light radiance must be rejected.");

    request = PrimaryRenderRequest{};
    request.directionalLight.incidentRadiance.z =
        std::numeric_limits<float>::infinity();
    expectInvalidArgument(
        [&] { request.validate(); }, "Non-finite directional-light radiance must be rejected.");

    request = PrimaryRenderRequest{};
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

void testDirectDiffuse() {
    const std::filesystem::path sourceDirectory{RAYM0NADE_TEST_SOURCE_DIR};
    Model model{sourceDirectory / "tests" / "data", "directional_light.obj", "null"};
    expect(model.faceCount() == 2U, "The directional-light fixture must contain two triangles.");

    PrimaryRenderRequest request{};
    request.extent = ImageExtent{13U, 9U};
    request.camera.pixelScale = 0.10F;
    request.directionalLight.directionToLight = vec3{3.0F, 0.0F, -4.0F};
    request.directionalLight.incidentRadiance = vec3{kPi, 2.0F * kPi, 0.5F * kPi};
    request.aov = PrimaryAov::DirectDiffuse;

    const LinearImage direct = renderPrimaryAovCpu(model, request);
    direct.validate();
    expectNear(
        pixelAt(direct, 0U, 0U),
        vec3{0.0F},
        0.0F,
        "A direct-diffuse miss must remain black.");
    expectNear(
        pixelAt(direct, 4U, 4U),
        vec3{0.64F, 0.40F, 0.04F},
        1.0e-6F,
        "A visible back-facing receiver must use the oriented normal and Lambert cosine.");
    expectNear(
        pixelAt(direct, 6U, 4U),
        vec3{0.0F},
        0.0F,
        "The blocker must cast a deterministic hard shadow on the receiver.");
    expectNear(
        pixelAt(direct, 9U, 4U),
        vec3{0.12F, 0.80F, 0.36F},
        1.0e-6F,
        "The visible blocker must receive the same normalized directional light.");

    const LinearImage repeated = renderPrimaryAovCpu(model, request);
    expectImagesExactlyEqual(
        repeated, direct, "Repeated CPU direct-diffuse renders must be pixel-identical.");

    PrimaryRenderRequest scaledDirection = request;
    scaledDirection.directionalLight.directionToLight = vec3{0.75F, 0.0F, -1.0F};
    const LinearImage normalized = renderPrimaryAovCpu(model, scaledDirection);
    expectImagesExactlyEqual(
        normalized,
        direct,
        "Scaling a directional-light vector must not change the direct-diffuse image.");

    PrimaryRenderRequest backLight = request;
    backLight.directionalLight.directionToLight = vec3{0.0F, 0.0F, 1.0F};
    const LinearImage unlit = renderPrimaryAovCpu(model, backLight);
    for (const vec3& pixel : unlit.pixels) {
        expect(
            exactlyEqual(pixel, vec3{0.0F}),
            "A light behind every camera-facing shape normal must contribute nothing.");
    }

    PrimaryRenderRequest twoSided{};
    twoSided.extent = ImageExtent{4U, 4U};
    twoSided.camera.pixelScale = 0.20F;
    twoSided.directionalLight.directionToLight = vec3{0.0F, 0.0F, -1.0F};
    twoSided.directionalLight.incidentRadiance = vec3{kPi};
    twoSided.aov = PrimaryAov::DirectDiffuse;
    const LinearImage viewedFromBelow = renderPrimaryAovCpu(model, twoSided);
    expectNear(
        pixelAt(viewedFromBelow, 2U, 2U),
        vec3{0.8F, 0.25F, 0.1F},
        1.0e-6F,
        "The receiver back face must be lit after orienting its shape normal toward the camera.");

    twoSided.camera.position = vec3{0.0F, 0.0F, 6.0F};
    twoSided.camera.direction = vec3{0.0F, 0.0F, -1.0F};
    twoSided.directionalLight.directionToLight = vec3{0.0F, 0.0F, 1.0F};
    const LinearImage viewedFromAbove = renderPrimaryAovCpu(model, twoSided);
    expectNear(
        pixelAt(viewedFromAbove, 2U, 2U),
        vec3{0.8F, 0.25F, 0.1F},
        1.0e-6F,
        "The receiver front face must retain the same two-sided direct-diffuse response.");
}

void testCpuPrimaryThreading() {
    const std::filesystem::path sourceDirectory{RAYM0NADE_TEST_SOURCE_DIR};
    Model model{
        sourceDirectory / "tests" / "data", "directional_light.obj", "null"};

    PrimaryRenderRequest request{};
    request.extent = ImageExtent{13U, 9U};
    request.camera.pixelScale = 0.10F;
    request.directionalLight.directionToLight = vec3{3.0F, 0.0F, -4.0F};
    request.directionalLight.incidentRadiance =
        vec3{kPi, 2.0F * kPi, 0.5F * kPi};

    constexpr std::array<PrimaryAov, 3> aovs{
        PrimaryAov::BaseColor,
        PrimaryAov::ShapeNormal,
        PrimaryAov::DirectDiffuse,
    };
    constexpr std::array<CpuPrimaryRenderOptions, 4> options{{
        CpuPrimaryRenderOptions{1},
        CpuPrimaryRenderOptions{2},
        CpuPrimaryRenderOptions{4},
        CpuPrimaryRenderOptions{0},
    }};
    expect(
        CpuPrimaryRenderOptions{}.resolvedThreadCount(9U) == 1U,
        "Default CPU primary-render options must preserve the single-thread oracle.");
    expect(
        CpuPrimaryRenderOptions{2}.resolvedThreadCount(9U) == 2U &&
            CpuPrimaryRenderOptions{4}.resolvedThreadCount(9U) == 4U,
        "Explicit CPU primary-render worker counts must be preserved.");
    expect(
        CpuPrimaryRenderOptions{4}.resolvedThreadCount(1U) == 1U,
        "CPU primary-render worker count must be clamped to the image height.");
    const std::uint32_t automaticWorkers =
        CpuPrimaryRenderOptions{0}.resolvedThreadCount(9U);
    expect(
        automaticWorkers >= 1U && automaticWorkers <= 9U,
        "Automatic CPU primary-render worker selection must be valid.");

    for (PrimaryAov aov : aovs) {
        request.aov = aov;
        const LinearImage singleThreadReference =
            renderPrimaryAovCpu(model, request);
        for (const CpuPrimaryRenderOptions& option : options) {
            const LinearImage threaded =
                renderPrimaryAovCpu(model, request, option);
            expectImagesExactlyEqual(
                threaded,
                singleThreadReference,
                "CPU primary AOV output changed with thread count " +
                    std::to_string(option.threadCount) + '.');
        }
    }

    CpuPrimaryRenderOptions invalidOptions;
    invalidOptions.threadCount = -1;
    expectInvalidArgument(
        [&] { (void)renderPrimaryAovCpu(model, request, invalidOptions); },
        "A negative CPU primary-render thread count must be rejected.");
}

}  // namespace

int main() {
    try {
        testContractValidation();
        testLegacyRayConvention();
        testPrimaryAovs();
        testDirectDiffuse();
        testCpuPrimaryThreading();
        std::cout << "Primary render contract test passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
