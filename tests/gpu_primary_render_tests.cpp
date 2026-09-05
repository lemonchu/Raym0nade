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
constexpr float kTexturePixelTolerance = 2.0e-5F;
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

constexpr std::uint32_t packRgba8(
    std::uint32_t red,
    std::uint32_t green,
    std::uint32_t blue,
    std::uint32_t alpha) noexcept {
    return red | (green << 8U) | (blue << 16U) | (alpha << 24U);
}

PackedVertex makePackedVertex(const vec3& position, const vec2& uv) noexcept {
    return PackedVertex{
        {position.x, position.y, position.z, 0.0F},
        {0.0F, -1.0F, uv.x, uv.y},
    };
}

PackedMaterial makePackedMaterial(const vec3& diffuse) noexcept {
    PackedMaterial material;
    material.diffuseAndOpacity = {diffuse.x, diffuse.y, diffuse.z, 1.0F};
    material.emissionAndIor = {0.0F, 0.0F, 0.0F, 1.5F};
    material.transmissionAndRoughness = {1.0F, 1.0F, 1.0F, 0.5F};
    material.metallicSpecularAndReserved = {0.0F, 0.04F, 0.0F, 0.0F};
    return material;
}

void appendTriangle(
    PackedSceneData& scene,
    const std::array<vec3, 3>& positions,
    const std::array<vec2, 3>& uvs,
    std::uint32_t materialId) {
    const std::uint32_t firstVertex =
        static_cast<std::uint32_t>(scene.vertices.size());
    for (std::size_t corner = 0U; corner < positions.size(); ++corner) {
        scene.vertices.push_back(makePackedVertex(positions[corner], uvs[corner]));
    }
    scene.triangleIndices.insert(
        scene.triangleIndices.end(),
        {firstVertex, firstVertex + 1U, firstVertex + 2U});
    scene.triangleMaterialIds.push_back(materialId);
}

void bindDiffuseTexture(
    PackedSceneData& scene,
    std::uint32_t materialId,
    std::uint32_t textureId,
    bool cutout) {
    PackedMaterial& material = scene.materials.at(materialId);
    material.flagsAndReserved[0] |= kPackedMaterialHasDiffuseTexture;
    if (cutout) {
        material.flagsAndReserved[0] |= kPackedMaterialCutout;
    }
    material.textureIds[0] = textureId;
}

void addEncodedFilterTexture(PackedSceneData& scene, std::uint32_t materialId) {
    const std::uint32_t textureId =
        static_cast<std::uint32_t>(scene.textures.size());
    const std::uint32_t firstMip =
        static_cast<std::uint32_t>(scene.textureMipLevels.size());
    const std::uint32_t firstTexel =
        static_cast<std::uint32_t>(scene.textureTexelsRgba8.size());
    scene.textures.push_back(PackedTexture{firstMip, 3U, 1U, 4U});
    scene.textureMipLevels.insert(
        scene.textureMipLevels.end(),
        {
            PackedTextureMip{firstTexel, 4U, 1U, 4U},
            PackedTextureMip{firstTexel + 4U, 2U, 1U, 2U},
            PackedTextureMip{firstTexel + 6U, 1U, 1U, 1U},
        });
    scene.textureTexelsRgba8.insert(
        scene.textureTexelsRgba8.end(),
        {
            packRgba8(0U, 0U, 0U, 255U),
            packRgba8(64U, 64U, 64U, 255U),
            packRgba8(192U, 192U, 192U, 255U),
            packRgba8(255U, 255U, 255U, 255U),
            packRgba8(32U, 32U, 32U, 255U),
            packRgba8(223U, 223U, 223U, 255U),
            packRgba8(127U, 127U, 127U, 255U),
        });
    bindDiffuseTexture(scene, materialId, textureId, false);
}

void addCutoutThresholdTexture(PackedSceneData& scene, std::uint32_t materialId) {
    const std::uint32_t textureId =
        static_cast<std::uint32_t>(scene.textures.size());
    const std::uint32_t firstMip =
        static_cast<std::uint32_t>(scene.textureMipLevels.size());
    const std::uint32_t firstTexel =
        static_cast<std::uint32_t>(scene.textureTexelsRgba8.size());
    scene.textures.push_back(PackedTexture{firstMip, 2U, 2U, 1U});
    scene.textureMipLevels.insert(
        scene.textureMipLevels.end(),
        {
            PackedTextureMip{firstTexel, 2U, 2U, 1U},
            PackedTextureMip{firstTexel + 2U, 1U, 1U, 1U},
        });
    scene.textureTexelsRgba8.insert(
        scene.textureTexelsRgba8.end(),
        {
            packRgba8(255U, 0U, 0U, 0U),
            packRgba8(255U, 0U, 0U, 1U),
            packRgba8(255U, 0U, 0U, 0U),
        });
    bindDiffuseTexture(scene, materialId, textureId, true);
}

vec3 encodedOrientedShapeNormal(
    const std::array<vec3, 3>& positions,
    const vec3& incomingDirection) noexcept {
    vec3 normal = safeNormalize(
        glm::cross(positions[1] - positions[0], positions[2] - positions[0]),
        -incomingDirection);
    if (glm::dot(normal, -incomingDirection) < 0.0F) {
        normal = -normal;
    }
    return (normal + vec3{1.0F}) * 0.5F;
}

void testTexturedBaseColor(
    const std::filesystem::path& shaderPath,
    const VulkanRayQueryOptions& options) {
    PackedSceneData scene;
    const vec3 diffuseFactor{0.8F, 0.4F, 0.2F};
    scene.materials.push_back(makePackedMaterial(diffuseFactor));
    appendTriangle(
        scene,
        {
            vec3{-2.0F, -2.0F, 2.0F},
            vec3{2.0F, -2.0F, 2.0F},
            vec3{0.0F, 2.0F, 2.0F},
        },
        {
            vec2{0.0F, 1.5F},
            vec2{0.0F, 2.0F},
            vec2{0.0F, 1.875F},
        },
        0U);
    addEncodedFilterTexture(scene, 0U);
    scene.validate();

    VulkanPrimaryRenderer renderer{scene, shaderPath, options};
    PrimaryRenderRequest request;
    request.extent = ImageExtent{2U, 2U};
    request.camera.pixelScale = 0.10F;
    request.aov = PrimaryAov::BaseColor;

    const VulkanPrimaryRenderResult first = renderer.render(request);
    const VulkanPrimaryRenderResult second = renderer.render(request);
    constexpr float encodedBilinearValue = 48.0F / 255.0F;
    const vec3 cpuExpected =
        vec3{std::pow(encodedBilinearValue, 2.2F)} * diffuseFactor;
    expectNear(
        pixelAt(first.image, 1U, 1U),
        cpuExpected,
        kTexturePixelTolerance,
        "Textured BaseColor must wrap and flip V, filter encoded texels, then apply gamma and the diffuse factor.");
    compareRepeatedGpuRender(first.image, second.image, "Textured BaseColor");
    validateTimings(first.timings, "Textured BaseColor first render");
    validateTimings(second.timings, "Textured BaseColor repeated render");
    expectCleanValidation(renderer.validationReport(), "Textured BaseColor renderer");
}

void testPrimaryCutout(
    const std::filesystem::path& shaderPath,
    const VulkanRayQueryOptions& options) {
    PackedSceneData scene;
    const vec3 backDiffuse{0.1F, 0.7F, 0.3F};
    scene.materials.push_back(makePackedMaterial(vec3{1.0F}));
    scene.materials.push_back(makePackedMaterial(backDiffuse));

    const std::array<vec3, 3> frontPositions{
        vec3{-2.0F, -2.0F, 2.0F},
        vec3{2.0F, -2.0F, 2.0F},
        vec3{0.0F, 2.0F, 2.0F},
    };
    appendTriangle(
        scene,
        frontPositions,
        {
            vec2{-0.75F, 0.0F},
            vec2{1.25F, 0.0F},
            vec2{0.25F, 0.0F},
        },
        0U);
    const std::array<vec3, 3> backPositions{
        vec3{-3.0F, -2.0F, 3.0F},
        vec3{3.0F, -2.0F, 3.0F},
        vec3{0.0F, 2.0F, 4.0F},
    };
    appendTriangle(
        scene,
        backPositions,
        {vec2{0.0F}, vec2{0.0F}, vec2{0.0F}},
        1U);
    addCutoutThresholdTexture(scene, 0U);
    scene.validate();

    VulkanPrimaryRenderer renderer{scene, shaderPath, options};
    PrimaryRenderRequest request;
    request.extent = ImageExtent{4U, 2U};
    request.camera.pixelScale = 0.25F;
    request.aov = PrimaryAov::BaseColor;

    const VulkanPrimaryRenderResult baseColorFirst = renderer.render(request);
    const VulkanPrimaryRenderResult baseColorSecond = renderer.render(request);
    expectNear(
        pixelAt(baseColorFirst.image, 1U, 1U),
        backDiffuse,
        kTexturePixelTolerance,
        "A zero-alpha primary candidate must be ignored in favor of the rear triangle.");
    expectNear(
        pixelAt(baseColorFirst.image, 3U, 1U),
        vec3{1.0F, 0.0F, 0.0F},
        kTexturePixelTolerance,
        "A diffuse alpha of 1/255 must exceed the cutout threshold and commit the front triangle.");
    compareRepeatedGpuRender(
        baseColorFirst.image, baseColorSecond.image, "Cutout BaseColor");

    request.aov = PrimaryAov::ShapeNormal;
    const VulkanPrimaryRenderResult shapeNormalFirst = renderer.render(request);
    const VulkanPrimaryRenderResult shapeNormalSecond = renderer.render(request);
    const vec3 expectedBackNormal = encodedOrientedShapeNormal(
        backPositions, makePrimaryRay(request.camera, request.extent, 1U, 1U).direction);
    const vec3 expectedFrontNormal = encodedOrientedShapeNormal(
        frontPositions, makePrimaryRay(request.camera, request.extent, 3U, 1U).direction);
    expectNear(
        pixelAt(shapeNormalFirst.image, 1U, 1U),
        expectedBackNormal,
        kTexturePixelTolerance,
        "ShapeNormal must describe the rear triangle after rejecting a zero-alpha candidate.");
    expectNear(
        pixelAt(shapeNormalFirst.image, 3U, 1U),
        expectedFrontNormal,
        kTexturePixelTolerance,
        "ShapeNormal must describe the front triangle when diffuse alpha is 1/255.");
    compareRepeatedGpuRender(
        shapeNormalFirst.image, shapeNormalSecond.image, "Cutout ShapeNormal");

    validateTimings(baseColorFirst.timings, "Cutout BaseColor first render");
    validateTimings(baseColorSecond.timings, "Cutout BaseColor repeated render");
    validateTimings(shapeNormalFirst.timings, "Cutout ShapeNormal first render");
    validateTimings(shapeNormalSecond.timings, "Cutout ShapeNormal repeated render");
    expectCleanValidation(renderer.validationReport(), "Primary-cutout renderer");
}

void testDirectDiffuseCutoutShadow(
    const std::filesystem::path& shaderPath,
    const VulkanRayQueryOptions& options) {
    PackedSceneData scene;
    const vec3 receiverDiffuse{0.25F, 0.5F, 0.75F};
    scene.materials.push_back(makePackedMaterial(receiverDiffuse));
    scene.materials.push_back(makePackedMaterial(vec3{1.0F}));
    appendTriangle(
        scene,
        {
            vec3{-4.0F, -3.0F, 3.0F},
            vec3{4.0F, -3.0F, 3.0F},
            vec3{0.0F, 3.0F, 3.0F},
        },
        {vec2{0.0F}, vec2{0.0F}, vec2{0.0F}},
        0U);
    appendTriangle(
        scene,
        {
            vec3{0.10F, -0.25F, 2.0F},
            vec3{0.40F, -0.25F, 2.0F},
            vec3{0.25F, 0.25F, 2.0F},
        },
        {vec2{0.0F}, vec2{0.0F}, vec2{0.0F}},
        1U);
    appendTriangle(
        scene,
        {
            vec3{1.60F, -0.25F, 2.0F},
            vec3{1.90F, -0.25F, 2.0F},
            vec3{1.75F, 0.25F, 2.0F},
        },
        {vec2{0.5F, 0.0F}, vec2{0.5F, 0.0F}, vec2{0.5F, 0.0F}},
        1U);
    addCutoutThresholdTexture(scene, 1U);
    scene.validate();

    VulkanPrimaryRenderer renderer{scene, shaderPath, options};
    PrimaryRenderRequest request;
    request.extent = ImageExtent{4U, 2U};
    request.camera.pixelScale = 0.25F;
    request.directionalLight.directionToLight = vec3{1.0F, 0.0F, -1.0F};
    const float radianceScale = kPi * std::sqrt(2.0F);
    request.directionalLight.incidentRadiance =
        vec3{radianceScale, 0.5F * radianceScale, 2.0F * radianceScale};
    request.aov = PrimaryAov::DirectDiffuse;

    const VulkanPrimaryRenderResult first = renderer.render(request);
    const VulkanPrimaryRenderResult second = renderer.render(request);
    const vec3 directionToLight =
        safeNormalize(request.directionalLight.directionToLight);
    const vec3 expectedLit = expectedDirectDiffuse(
        receiverDiffuse,
        request.directionalLight.incidentRadiance,
        glm::dot(vec3{0.0F, 0.0F, -1.0F}, directionToLight));
    expectNearAbsoluteRelative(
        pixelAt(first.image, 1U, 1U),
        expectedLit,
        kDirectDiffuseAbsoluteTolerance,
        kDirectDiffuseRelativeTolerance,
        "A zero-alpha blocker must not shadow the DirectDiffuse receiver.");
    expectExactlyZeroPixel(
        first.image,
        3U,
        1U,
        "A blocker with diffuse alpha 1/255 must cast an exact-zero DirectDiffuse shadow.");
    compareRepeatedGpuRender(first.image, second.image, "Cutout-shadow DirectDiffuse");
    validateTimings(first.timings, "Cutout-shadow DirectDiffuse first render");
    validateTimings(second.timings, "Cutout-shadow DirectDiffuse repeated render");
    expectCleanValidation(renderer.validationReport(), "Cutout-shadow renderer");
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

    PackedSceneData halfOpacityScene = scene;
    for (std::uint32_t materialId : halfOpacityScene.triangleMaterialIds) {
        halfOpacityScene.materials[materialId].diffuseAndOpacity[3] = 0.5F;
    }
    {
        VulkanPrimaryRenderer halfOpacityRenderer{
            halfOpacityScene, shaderPath, options};
        const VulkanPrimaryRenderResult halfOpacityBaseColor =
            halfOpacityRenderer.render(baseColorRequest);
        const VulkanPrimaryRenderResult halfOpacityDirect =
            halfOpacityRenderer.render(directRequest);
        compareWithCpu(
            cpuBaseColor,
            halfOpacityBaseColor.image,
            "Half-opacity BaseColor");
        compareDirectDiffuseWithCpu(
            cpuDirect,
            halfOpacityDirect.image,
            "Half-opacity DirectDiffuse");
        validateTimings(
            halfOpacityBaseColor.timings, "Half-opacity BaseColor render");
        validateTimings(
            halfOpacityDirect.timings, "Half-opacity DirectDiffuse render");
        expectCleanValidation(
            halfOpacityRenderer.validationReport(), "Half-opacity renderer");
    }

    PackedSceneData zeroOpacityScene = scene;
    const vec3 transmissionColor{0.35F, -0.25F, 0.55F};
    for (std::uint32_t materialId : zeroOpacityScene.triangleMaterialIds) {
        PackedMaterial& material = zeroOpacityScene.materials[materialId];
        material.diffuseAndOpacity[3] = 0.0F;
        material.transmissionAndRoughness[0] = transmissionColor.x;
        material.transmissionAndRoughness[1] = transmissionColor.y;
        material.transmissionAndRoughness[2] = transmissionColor.z;
    }
    {
        VulkanPrimaryRenderer zeroOpacityRenderer{
            zeroOpacityScene, shaderPath, options};
        const VulkanPrimaryRenderResult zeroOpacityBaseColor =
            zeroOpacityRenderer.render(baseColorRequest);
        for (std::size_t index = 0U; index < cpuBaseColor.pixels.size(); ++index) {
            const vec3 expected = exactlyEqual(cpuBaseColor.pixels[index], vec3{0.0F})
                                      ? vec3{0.0F}
                                      : transmissionColor;
            expectNear(
                zeroOpacityBaseColor.image.pixels[index],
                expected,
                kPixelTolerance,
                "Zero-opacity BaseColor must use the scalar transmission color.");
        }

        const VulkanPrimaryRenderResult zeroOpacityDirect =
            zeroOpacityRenderer.render(directRequest);
        const vec3 clampedTransmission{transmissionColor.x, 0.0F, transmissionColor.z};
        const vec3 expectedTransmissionDirect = expectedDirectDiffuse(
            clampedTransmission,
            directRequest.directionalLight.incidentRadiance,
            0.8F);
        expectNearAbsoluteRelative(
            pixelAt(zeroOpacityDirect.image, 4U, 4U),
            expectedTransmissionDirect,
            kDirectDiffuseAbsoluteTolerance,
            kDirectDiffuseRelativeTolerance,
            "Zero-opacity DirectDiffuse must clamp transmission before Lambert shading.");
        expectNearAbsoluteRelative(
            pixelAt(zeroOpacityDirect.image, 9U, 4U),
            expectedTransmissionDirect,
            kDirectDiffuseAbsoluteTolerance,
            kDirectDiffuseRelativeTolerance,
            "A visible zero-opacity blocker must use clamped transmission for DirectDiffuse.");
        expectExactlyZeroPixel(
            zeroOpacityDirect.image,
            6U,
            4U,
            "Scalar zero opacity must not disable hard-shadow occlusion.");
        validateTimings(
            zeroOpacityBaseColor.timings, "Zero-opacity BaseColor render");
        validateTimings(
            zeroOpacityDirect.timings, "Zero-opacity DirectDiffuse render");
        expectCleanValidation(
            zeroOpacityRenderer.validationReport(), "Zero-opacity renderer");
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
        const VulkanPrimaryRenderResult texturedFirst =
            texturedRenderer.render(directRequest);
        const VulkanPrimaryRenderResult texturedSecond =
            texturedRenderer.render(directRequest);
        compareDirectDiffuseWithCpu(
            cpuDirect,
            texturedFirst.image,
            "White-textured DirectDiffuse");
        compareDirectDiffuseWithCpu(
            cpuDirect,
            texturedSecond.image,
            "Repeated white-textured DirectDiffuse");
        compareRepeatedGpuRender(
            texturedFirst.image,
            texturedSecond.image,
            "White-textured DirectDiffuse");
        validateTimings(
            texturedFirst.timings, "White-textured DirectDiffuse first render");
        validateTimings(
            texturedSecond.timings, "White-textured DirectDiffuse repeated render");
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

    testTexturedBaseColor(shaderPath, options);
    testPrimaryCutout(shaderPath, options);
    testDirectDiffuseCutoutShadow(shaderPath, options);
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
