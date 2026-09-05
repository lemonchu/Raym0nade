#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "raym0nade/gpu/vulkan_primary_renderer.hpp"
#include "raym0nade/model.hpp"
#include "raym0nade/render.hpp"
#include "raym0nade/scene_data.hpp"

#ifndef RAYM0NADE_TEST_SOURCE_DIR
#error "RAYM0NADE_TEST_SOURCE_DIR must be defined by CMake."
#endif

#ifndef RAYM0NADE_PRIMARY_AOV_SHADER
#error "RAYM0NADE_PRIMARY_AOV_SHADER must be defined by CMake."
#endif

namespace {

using namespace raym0nade;
using namespace raym0nade::gpu;

constexpr float kPixelTolerance = 1.0e-6F;
constexpr std::string_view kNoCompatibleDevice =
    "No AMD Vulkan device satisfies the Ray Query backend requirements.";
constexpr std::array<const char*, 4> kExpectedHitMask{
    "1111",
    "0111",
    "0111",
    "0010",
};

void expect(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool exactlyEqual(const vec3& left, const vec3& right) noexcept {
    return left.x == right.x && left.y == right.y && left.z == right.z;
}

void expectNear(
    const vec3& actual,
    const vec3& expected,
    float tolerance,
    const std::string& message) {
    if (!isFinite(actual) ||
        std::abs(actual.x - expected.x) > tolerance ||
        std::abs(actual.y - expected.y) > tolerance ||
        std::abs(actual.z - expected.z) > tolerance) {
        throw std::runtime_error(message);
    }
}

void validateReferenceFixture(const LinearImage& image, PrimaryAov aov) {
    image.validate();
    for (std::uint32_t y = 0U; y < image.extent.height; ++y) {
        for (std::uint32_t x = 0U; x < image.extent.width; ++x) {
            const std::size_t index =
                static_cast<std::size_t>(y) * image.extent.width + x;
            const bool expectedHit = kExpectedHitMask[y][x] == '1';
            const vec3 expected =
                !expectedHit
                    ? vec3{0.0F}
                    : (aov == PrimaryAov::BaseColor
                           ? vec3{0.8F, 0.2F, 0.1F}
                           : vec3{0.5F, 0.5F, 0.0F});
            expectNear(
                image.pixels[index],
                expected,
                kPixelTolerance,
                "The CPU reference fixture produced an unexpected primary AOV pixel.");
        }
    }
}

void compareWithCpu(
    const LinearImage& cpu,
    const LinearImage& gpu,
    const std::string& aovName) {
    cpu.validate();
    gpu.validate();
    expect(
        cpu.extent.width == gpu.extent.width &&
            cpu.extent.height == gpu.extent.height,
        aovName + " CPU/GPU extents differ.");
    expect(cpu.pixels.size() == gpu.pixels.size(), aovName + " CPU/GPU sizes differ.");
    for (std::size_t index = 0U; index < cpu.pixels.size(); ++index) {
        expectNear(
            gpu.pixels[index],
            cpu.pixels[index],
            kPixelTolerance,
            aovName + " CPU/GPU pixels differ at index " + std::to_string(index) + '.');
    }
}

void compareRepeatedGpuRender(
    const LinearImage& first,
    const LinearImage& second,
    const std::string& aovName) {
    first.validate();
    second.validate();
    expect(
        first.extent.width == second.extent.width &&
            first.extent.height == second.extent.height &&
            first.pixels.size() == second.pixels.size(),
        aovName + " repeated GPU render dimensions changed.");
    for (std::size_t index = 0U; index < first.pixels.size(); ++index) {
        expect(
            exactlyEqual(first.pixels[index], second.pixels[index]),
            aovName + " repeated GPU renders differ at index " +
                std::to_string(index) + '.');
    }
}

void validateTimings(
    const VulkanPrimaryRenderTimings& timings, const std::string& description) {
    expect(
        std::isfinite(timings.dispatchAndReadbackMilliseconds) &&
            timings.dispatchAndReadbackMilliseconds >= 0.0,
        description + " host dispatch/readback timing is invalid.");
    if (timings.gpuTimestampAvailable) {
        expect(
            std::isfinite(timings.gpuDispatchMilliseconds) &&
                timings.gpuDispatchMilliseconds >= 0.0,
            description + " GPU dispatch timing is invalid.");
    }
}

void printTimings(
    const std::string& description, const VulkanPrimaryRenderTimings& timings) {
    std::cout << description << ": host dispatch/readback "
              << timings.dispatchAndReadbackMilliseconds << " ms, GPU dispatch ";
    if (timings.gpuTimestampAvailable) {
        std::cout << timings.gpuDispatchMilliseconds << " ms\n";
    } else {
        std::cout << "unavailable\n";
    }
}

int runTest() {
    const std::filesystem::path sourceDirectory{RAYM0NADE_TEST_SOURCE_DIR};
    Model model{sourceDirectory / "tests" / "data", "triangle.obj", "null"};
    const PackedSceneData scene = model.packScene();

    VulkanRayQueryOptions options;
    options.requestValidation = true;
    VulkanPrimaryRenderer renderer{
        scene, std::filesystem::path{RAYM0NADE_PRIMARY_AOV_SHADER}, options};

    PrimaryRenderRequest request;
    request.extent = ImageExtent{4U, 4U};
    request.camera.pixelScale = 0.20F;
    request.aov = PrimaryAov::BaseColor;

    const LinearImage cpuBaseColor = renderPrimaryAovCpu(model, request);
    validateReferenceFixture(cpuBaseColor, request.aov);
    const VulkanPrimaryRenderResult gpuBaseColorFirst = renderer.render(request);
    const VulkanPrimaryRenderResult gpuBaseColorSecond = renderer.render(request);
    compareWithCpu(cpuBaseColor, gpuBaseColorFirst.image, "BaseColor");
    compareWithCpu(cpuBaseColor, gpuBaseColorSecond.image, "BaseColor");
    compareRepeatedGpuRender(
        gpuBaseColorFirst.image, gpuBaseColorSecond.image, "BaseColor");
    validateTimings(gpuBaseColorFirst.timings, "BaseColor first render");
    validateTimings(gpuBaseColorSecond.timings, "BaseColor second render");

    request.aov = PrimaryAov::ShapeNormal;
    const LinearImage cpuShapeNormal = renderPrimaryAovCpu(model, request);
    validateReferenceFixture(cpuShapeNormal, request.aov);
    const VulkanPrimaryRenderResult gpuShapeNormalFirst = renderer.render(request);
    const VulkanPrimaryRenderResult gpuShapeNormalSecond = renderer.render(request);
    compareWithCpu(cpuShapeNormal, gpuShapeNormalFirst.image, "ShapeNormal");
    compareWithCpu(cpuShapeNormal, gpuShapeNormalSecond.image, "ShapeNormal");
    compareRepeatedGpuRender(
        gpuShapeNormalFirst.image, gpuShapeNormalSecond.image, "ShapeNormal");
    validateTimings(gpuShapeNormalFirst.timings, "ShapeNormal first render");
    validateTimings(gpuShapeNormalSecond.timings, "ShapeNormal second render");

    PrimaryRenderRequest expandedRequest;
    expandedRequest.extent = ImageExtent{13U, 9U};
    expandedRequest.camera.position = vec3{0.037F, -0.023F, 0.0F};
    expandedRequest.camera.right = vec3{1.30F, 0.0F, 0.0F};
    expandedRequest.camera.up = vec3{0.0F, 0.80F, 0.0F};
    expandedRequest.camera.pixelScale = 0.08F;
    expandedRequest.aov = PrimaryAov::BaseColor;

    const LinearImage cpuExpanded = renderPrimaryAovCpu(model, expandedRequest);
    std::size_t expandedHitCount = 0U;
    for (const vec3& pixel : cpuExpanded.pixels) {
        if (!exactlyEqual(pixel, vec3{0.0F})) {
            ++expandedHitCount;
        }
    }
    expect(
        expandedHitCount > 0U && expandedHitCount < cpuExpanded.pixels.size(),
        "The expanded CPU reference must contain stable hit and miss pixels.");
    const VulkanPrimaryRenderResult gpuExpanded = renderer.render(expandedRequest);
    compareWithCpu(cpuExpanded, gpuExpanded.image, "Expanded BaseColor");
    validateTimings(gpuExpanded.timings, "Expanded BaseColor render");

    request.aov = PrimaryAov::BaseColor;
    const VulkanPrimaryRenderResult gpuBaseColorAfterGrowth = renderer.render(request);
    compareWithCpu(cpuBaseColor, gpuBaseColorAfterGrowth.image, "BaseColor after growth");
    compareRepeatedGpuRender(
        gpuBaseColorFirst.image, gpuBaseColorAfterGrowth.image, "BaseColor after growth");
    validateTimings(
        gpuBaseColorAfterGrowth.timings, "BaseColor render after output-buffer growth");

    const VulkanValidationReport validation = renderer.validationReport();
    if (validation.errorCount != 0U || validation.warningCount != 0U) {
        throw std::runtime_error(
            "Vulkan validation reported " + std::to_string(validation.errorCount) +
            " error(s) and " + std::to_string(validation.warningCount) + " warning(s).");
    }

    const VulkanRayQuerySetupTimings& setup = renderer.setupTimings();
    std::cout << std::fixed << std::setprecision(3)
              << "GPU primary renderer device: " << renderer.deviceName() << '\n'
              << "Setup: upload " << setup.uploadMilliseconds << " ms, AS build "
              << setup.accelerationStructureBuildMilliseconds << " ms\n";
    printTimings("BaseColor first", gpuBaseColorFirst.timings);
    printTimings("BaseColor second", gpuBaseColorSecond.timings);
    printTimings("ShapeNormal first", gpuShapeNormalFirst.timings);
    printTimings("ShapeNormal second", gpuShapeNormalSecond.timings);
    printTimings("Expanded BaseColor", gpuExpanded.timings);
    printTimings("BaseColor after growth", gpuBaseColorAfterGrowth.timings);
    std::cout << "Validation: requested=" << (validation.requested ? "yes" : "no")
              << ", layer=" << (validation.enabled ? "enabled" : "unavailable")
              << ", synchronization="
              << (validation.synchronizationValidationEnabled ? "enabled" : "unavailable")
              << " (" << validation.errorCount << " error(s), "
              << validation.warningCount << " warning(s))\n"
              << "CPU/GPU primary AOV test passed.\n";
    for (const std::string& message : validation.messages) {
        std::cout << "Validation message: " << message << '\n';
    }
    return 0;
}

}  // namespace

int main() {
    try {
        return runTest();
    } catch (const std::exception& error) {
        const std::string_view message{error.what()};
        if (message.rfind(kNoCompatibleDevice, 0U) == 0U) {
            std::cout << "SKIPPED: " << message << '\n';
            return 77;
        }
        std::cerr << "FAILED: " << message << '\n';
        return 1;
    }
}
