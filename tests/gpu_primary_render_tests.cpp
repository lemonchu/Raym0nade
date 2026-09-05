#include <algorithm>
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
constexpr float kDirectDiffuseAbsoluteTolerance = 1.0e-5F;
constexpr float kDirectDiffuseRelativeTolerance = 5.0e-5F;
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

template <typename Function>
void expectInvalidArgument(Function&& function, const std::string& message) {
    try {
        function();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

bool exactlyEqual(const vec3& left, const vec3& right) noexcept {
    return left.x == right.x && left.y == right.y && left.z == right.z;
}

const vec3& pixelAt(const LinearImage& image, std::uint32_t x, std::uint32_t y) {
    expect(
        x < image.extent.width && y < image.extent.height,
        "Test pixel coordinate is outside the image extent.");
    return image.pixels[static_cast<std::size_t>(y) * image.extent.width + x];
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

void expectNearAbsoluteRelative(
    const vec3& actual,
    const vec3& expected,
    float absoluteTolerance,
    float relativeTolerance,
    const std::string& message) {
    if (!isFinite(actual) || !isFinite(expected)) {
        throw std::runtime_error(message);
    }
    const auto componentMatches = [absoluteTolerance, relativeTolerance](
                                      float actualComponent, float expectedComponent) {
        const float scale = std::max(
            std::abs(actualComponent), std::abs(expectedComponent));
        return std::abs(actualComponent - expectedComponent) <=
               absoluteTolerance + relativeTolerance * scale;
    };
    if (!componentMatches(actual.x, expected.x) ||
        !componentMatches(actual.y, expected.y) ||
        !componentMatches(actual.z, expected.z)) {
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

void compareDirectDiffuseWithCpu(
    const LinearImage& cpu,
    const LinearImage& gpu,
    const std::string& description) {
    cpu.validate();
    gpu.validate();
    expect(
        cpu.extent.width == gpu.extent.width &&
            cpu.extent.height == gpu.extent.height,
        description + " CPU/GPU extents differ.");
    expect(
        cpu.pixels.size() == gpu.pixels.size(),
        description + " CPU/GPU sizes differ.");
    for (std::size_t index = 0U; index < cpu.pixels.size(); ++index) {
        expectNearAbsoluteRelative(
            gpu.pixels[index],
            cpu.pixels[index],
            kDirectDiffuseAbsoluteTolerance,
            kDirectDiffuseRelativeTolerance,
            description + " CPU/GPU pixels differ at index " +
                std::to_string(index) + '.');
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

void expectCleanValidation(
    const VulkanValidationReport& validation, const std::string& description) {
    if (validation.errorCount != 0U || validation.warningCount != 0U) {
        throw std::runtime_error(
            description + " reported " + std::to_string(validation.errorCount) +
            " Vulkan validation error(s) and " +
            std::to_string(validation.warningCount) + " warning(s).");
    }
}

void expectExactlyZeroPixel(
    const LinearImage& image,
    std::uint32_t x,
    std::uint32_t y,
    const std::string& message) {
    expect(exactlyEqual(pixelAt(image, x, y), vec3{0.0F}), message);
}

vec3 expectedDirectDiffuse(
    const vec3& diffuse,
    const vec3& incidentRadiance,
    float normalDotLight) noexcept {
    return diffuse * incidentRadiance * (normalDotLight / kPi);
}

void testDirectionalDirectDiffuse(
    const std::filesystem::path& sourceDirectory,
    const std::filesystem::path& shaderPath,
    const VulkanRayQueryOptions& options) {
    Model model{
        sourceDirectory / "tests" / "data", "directional_light.obj", "null"};
    expect(model.faceCount() == 2U, "The directional-light fixture must contain two triangles.");
    const PackedSceneData scene = model.packScene();
    VulkanPrimaryRenderer renderer{scene, shaderPath, options};

    const vec3 receiverDiffuse{0.8F, 0.25F, 0.1F};
    const vec3 blockerDiffuse{0.15F, 0.5F, 0.9F};

    PrimaryRenderRequest viewedFromBelow;
    viewedFromBelow.extent = ImageExtent{4U, 4U};
    viewedFromBelow.camera.pixelScale = 0.20F;
    viewedFromBelow.directionalLight.directionToLight = vec3{0.0F, 0.0F, -1.0F};
    viewedFromBelow.directionalLight.incidentRadiance = vec3{kPi};
    viewedFromBelow.aov = PrimaryAov::DirectDiffuse;
    const LinearImage cpuViewedFromBelow = renderPrimaryAovCpu(model, viewedFromBelow);
    const VulkanPrimaryRenderResult gpuViewedFromBelow = renderer.render(viewedFromBelow);
    compareDirectDiffuseWithCpu(
        cpuViewedFromBelow, gpuViewedFromBelow.image, "DirectDiffuse viewed from below");
    expectNear(
        pixelAt(cpuViewedFromBelow, 2U, 2U),
        receiverDiffuse,
        kPixelTolerance,
        "The CPU reference must light the receiver's camera-facing back side.");
    expectNear(
        pixelAt(gpuViewedFromBelow.image, 2U, 2U),
        receiverDiffuse,
        kPixelTolerance,
        "The GPU must light the receiver's camera-facing back side.");
    validateTimings(gpuViewedFromBelow.timings, "DirectDiffuse initial 4x4 render");

    PrimaryRenderRequest directRequest;
    directRequest.extent = ImageExtent{13U, 9U};
    directRequest.camera.pixelScale = 0.10F;
    directRequest.directionalLight.directionToLight = vec3{3.0F, 0.0F, -4.0F};
    directRequest.directionalLight.incidentRadiance =
        vec3{kPi, 2.0F * kPi, 0.5F * kPi};
    directRequest.aov = PrimaryAov::DirectDiffuse;

    const LinearImage cpuDirect = renderPrimaryAovCpu(model, directRequest);
    const VulkanPrimaryRenderResult gpuDirectFirst = renderer.render(directRequest);
    const VulkanPrimaryRenderResult gpuDirectSecond = renderer.render(directRequest);
    compareDirectDiffuseWithCpu(cpuDirect, gpuDirectFirst.image, "DirectDiffuse");
    compareDirectDiffuseWithCpu(
        cpuDirect, gpuDirectSecond.image, "Repeated DirectDiffuse");
    compareRepeatedGpuRender(
        gpuDirectFirst.image, gpuDirectSecond.image, "DirectDiffuse");

    const vec3 expectedReceiver = expectedDirectDiffuse(
        receiverDiffuse, directRequest.directionalLight.incidentRadiance, 0.8F);
    const vec3 expectedBlocker = expectedDirectDiffuse(
        blockerDiffuse, directRequest.directionalLight.incidentRadiance, 0.8F);
    expectNear(
        pixelAt(cpuDirect, 4U, 4U),
        expectedReceiver,
        kPixelTolerance,
        "The CPU clear receiver probe does not match the Lambert cosine contract.");
    expectNear(
        pixelAt(gpuDirectFirst.image, 4U, 4U),
        expectedReceiver,
        kPixelTolerance,
        "The GPU clear receiver probe does not match the Lambert cosine contract.");
    expectNear(
        pixelAt(cpuDirect, 9U, 4U),
        expectedBlocker,
        kPixelTolerance,
        "The CPU visible-blocker probe does not match the Lambert cosine contract.");
    expectNear(
        pixelAt(gpuDirectFirst.image, 9U, 4U),
        expectedBlocker,
        kPixelTolerance,
        "The GPU visible-blocker probe does not match the Lambert cosine contract.");
    expectExactlyZeroPixel(cpuDirect, 0U, 0U, "A CPU direct-diffuse miss must be exact zero.");
    expectExactlyZeroPixel(
        gpuDirectFirst.image, 0U, 0U, "A GPU direct-diffuse miss must be exact zero.");
    expectExactlyZeroPixel(
        cpuDirect, 6U, 4U, "The CPU blocker must produce an exact-zero hard shadow.");
    expectExactlyZeroPixel(
        gpuDirectFirst.image,
        6U,
        4U,
        "The GPU blocker must produce an exact-zero hard shadow.");
    validateTimings(gpuDirectFirst.timings, "DirectDiffuse first 13x9 render");
    validateTimings(gpuDirectSecond.timings, "DirectDiffuse repeated 13x9 render");

    PrimaryRenderRequest scaledDirection = directRequest;
    scaledDirection.directionalLight.directionToLight = vec3{0.75F, 0.0F, -1.0F};
    const LinearImage cpuScaledDirection = renderPrimaryAovCpu(model, scaledDirection);
    const VulkanPrimaryRenderResult gpuScaledDirection = renderer.render(scaledDirection);
    compareDirectDiffuseWithCpu(
        cpuScaledDirection, gpuScaledDirection.image, "Scaled-direction DirectDiffuse");
    compareRepeatedGpuRender(
        gpuDirectFirst.image,
        gpuScaledDirection.image,
        "DirectDiffuse direction-scale invariance");
    validateTimings(gpuScaledDirection.timings, "Scaled-direction DirectDiffuse render");

    PrimaryRenderRequest backLight = directRequest;
    backLight.directionalLight.directionToLight = vec3{0.0F, 0.0F, 1.0F};
    backLight.directionalLight.incidentRadiance = vec3{7.0F, 5.0F, 3.0F};
    const LinearImage cpuBackLight = renderPrimaryAovCpu(model, backLight);
    const VulkanPrimaryRenderResult gpuBackLight = renderer.render(backLight);
    compareDirectDiffuseWithCpu(
        cpuBackLight, gpuBackLight.image, "Back-light DirectDiffuse");
    for (std::uint32_t y = 0U; y < backLight.extent.height; ++y) {
        for (std::uint32_t x = 0U; x < backLight.extent.width; ++x) {
            expectExactlyZeroPixel(
                cpuBackLight,
                x,
                y,
                "A CPU light behind the oriented shape normal must contribute exact zero.");
            expectExactlyZeroPixel(
                gpuBackLight.image,
                x,
                y,
                "A GPU light behind the oriented shape normal must contribute exact zero.");
        }
    }
    validateTimings(gpuBackLight.timings, "Back-light DirectDiffuse render");

    PrimaryRenderRequest baseColorRequest = directRequest;
    baseColorRequest.directionalLight.directionToLight = vec3{-2.0F, 3.0F, 4.0F};
    baseColorRequest.directionalLight.incidentRadiance = vec3{9.0F, 8.0F, 7.0F};
    baseColorRequest.aov = PrimaryAov::BaseColor;
    const LinearImage cpuBaseColor = renderPrimaryAovCpu(model, baseColorRequest);
    const VulkanPrimaryRenderResult gpuBaseColor = renderer.render(baseColorRequest);
    compareWithCpu(cpuBaseColor, gpuBaseColor.image, "Directional fixture BaseColor");
    expectExactlyZeroPixel(
        cpuBaseColor, 0U, 0U, "The directional fixture miss classification changed.");
    expectNear(
        pixelAt(cpuBaseColor, 4U, 4U),
        receiverDiffuse,
        kPixelTolerance,
        "The clear probe must hit the receiver material.");
    expectNear(
        pixelAt(cpuBaseColor, 6U, 4U),
        receiverDiffuse,
        kPixelTolerance,
        "The shadow probe must first hit the receiver material.");
    expectNear(
        pixelAt(cpuBaseColor, 9U, 4U),
        blockerDiffuse,
        kPixelTolerance,
        "The blocker probe must first hit the blocker material.");
    validateTimings(gpuBaseColor.timings, "Directional fixture BaseColor render");

    PrimaryRenderRequest viewedFromAbove = viewedFromBelow;
    viewedFromAbove.camera.position = vec3{0.0F, 0.0F, 6.0F};
    viewedFromAbove.camera.direction = vec3{0.0F, 0.0F, -1.0F};
    viewedFromAbove.directionalLight.directionToLight = vec3{0.0F, 0.0F, 1.0F};
    const LinearImage cpuViewedFromAbove = renderPrimaryAovCpu(model, viewedFromAbove);
    const VulkanPrimaryRenderResult gpuViewedFromAbove = renderer.render(viewedFromAbove);
    compareDirectDiffuseWithCpu(
        cpuViewedFromAbove, gpuViewedFromAbove.image, "DirectDiffuse viewed from above");
    expectNear(
        pixelAt(gpuViewedFromAbove.image, 2U, 2U),
        receiverDiffuse,
        kPixelTolerance,
        "The GPU must preserve the receiver's front-side direct-diffuse response.");
    validateTimings(gpuViewedFromAbove.timings, "DirectDiffuse 4x4 render after growth");

    const VulkanPrimaryRenderResult gpuViewedFromBelowAfterChanges =
        renderer.render(viewedFromBelow);
    compareDirectDiffuseWithCpu(
        cpuViewedFromBelow,
        gpuViewedFromBelowAfterChanges.image,
        "DirectDiffuse viewed from below after push-constant changes");
    compareRepeatedGpuRender(
        gpuViewedFromBelow.image,
        gpuViewedFromBelowAfterChanges.image,
        "DirectDiffuse after output growth, shrink, and push-constant changes");
    validateTimings(
        gpuViewedFromBelowAfterChanges.timings,
        "DirectDiffuse repeated 4x4 render after growth");

    const VulkanValidationReport validation = renderer.validationReport();
    expectCleanValidation(validation, "Directional-light renderer");

    PackedSceneData transparentScene = scene;
    for (std::uint32_t materialId : transparentScene.triangleMaterialIds) {
        transparentScene.materials[materialId].diffuseAndOpacity[3] = 0.5F;
    }
    {
        VulkanPrimaryRenderer transparentRenderer{transparentScene, shaderPath, options};
        expectInvalidArgument(
            [&] { (void)transparentRenderer.render(directRequest); },
            "Vulkan DirectDiffuse must reject referenced transparent materials.");
        const VulkanValidationReport transparentValidation =
            transparentRenderer.validationReport();
        expectCleanValidation(
            transparentValidation, "Transparent-material capability test");
    }

    PackedSceneData texturedScene = scene;
    texturedScene.textures.push_back(PackedTexture{0U, 1U, 1U, 1U});
    texturedScene.textureMipLevels.push_back(
        PackedTextureMip{0U, 1U, 1U, 1U});
    texturedScene.textureTexelsRgba8.push_back(0xffffffffU);
    for (std::uint32_t materialId : texturedScene.triangleMaterialIds) {
        texturedScene.materials[materialId].flagsAndReserved[0] |=
            kPackedMaterialHasDiffuseTexture;
        texturedScene.materials[materialId].textureIds[0] = 0U;
    }
    {
        VulkanPrimaryRenderer texturedRenderer{texturedScene, shaderPath, options};
        expectInvalidArgument(
            [&] { (void)texturedRenderer.render(directRequest); },
            "Vulkan DirectDiffuse must reject referenced diffuse textures.");
        const VulkanValidationReport texturedValidation =
            texturedRenderer.validationReport();
        expectCleanValidation(
            texturedValidation, "Diffuse-texture capability test");
    }

    printTimings("DirectDiffuse initial 4x4", gpuViewedFromBelow.timings);
    printTimings("DirectDiffuse first 13x9", gpuDirectFirst.timings);
    printTimings("DirectDiffuse repeated 13x9", gpuDirectSecond.timings);
    printTimings("DirectDiffuse scaled direction", gpuScaledDirection.timings);
    printTimings("DirectDiffuse back light", gpuBackLight.timings);
    printTimings("Directional fixture BaseColor", gpuBaseColor.timings);
    printTimings("DirectDiffuse viewed from above", gpuViewedFromAbove.timings);
    printTimings(
        "DirectDiffuse final 4x4", gpuViewedFromBelowAfterChanges.timings);
    std::cout << "Directional renderer validation: requested="
              << (validation.requested ? "yes" : "no")
              << ", layer=" << (validation.enabled ? "enabled" : "unavailable")
              << ", synchronization="
              << (validation.synchronizationValidationEnabled ? "enabled" : "unavailable")
              << " (" << validation.errorCount << " error(s), "
              << validation.warningCount << " warning(s))\n";
}

int runTest() {
    const std::filesystem::path sourceDirectory{RAYM0NADE_TEST_SOURCE_DIR};
    const std::filesystem::path shaderPath{RAYM0NADE_PRIMARY_AOV_SHADER};
    Model model{sourceDirectory / "tests" / "data", "triangle.obj", "null"};
    const PackedSceneData scene = model.packScene();

    VulkanRayQueryOptions options;
    options.requestValidation = true;
    VulkanPrimaryRenderer renderer{scene, shaderPath, options};

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
    expectCleanValidation(validation, "Primary-AOV renderer");

    testDirectionalDirectDiffuse(sourceDirectory, shaderPath, options);

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
