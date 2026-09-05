#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "raym0nade/bvh.hpp"
#include "raym0nade/component.hpp"
#include "raym0nade/geometry.hpp"
#include "raym0nade/material.hpp"
#include "raym0nade/render.hpp"
#include "raym0nade/sampling.hpp"
#include "raym0nade/scene_data.hpp"

namespace {

using namespace raym0nade;

int failureCount = 0;

void expect(bool condition, std::string_view message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        ++failureCount;
    }
}

void expectNear(float actual, float expected, float tolerance, std::string_view message) {
    expect(std::abs(actual - expected) <= tolerance, message);
}

template <typename Function>
void expectInvalidArgument(Function&& function, std::string_view message) {
    try {
        std::forward<Function>(function)();
        expect(false, message);
    } catch (const std::invalid_argument&) {
    } catch (...) {
        expect(false, message);
    }
}

void expectUnitAndOrthogonal(
    const vec3& normal, const vec3& tangent, const vec3& bitangent, std::string_view context) {
    expect(isFinite(tangent) && isFinite(bitangent), context);
    expectNear(glm::length(tangent), 1.0F, 1.0e-5F, context);
    expectNear(glm::length(bitangent), 1.0F, 1.0e-5F, context);
    expectNear(glm::dot(normal, tangent), 0.0F, 1.0e-5F, context);
    expectNear(glm::dot(normal, bitangent), 0.0F, 1.0e-5F, context);
    expectNear(glm::dot(tangent, bitangent), 0.0F, 1.0e-5F, context);
}

Face triangleAt(float depth) {
    Face face;
    face.vertices[0] = vec3{-1.0F, -1.0F, depth};
    face.vertices[1] = vec3{1.0F, -1.0F, depth};
    face.vertices[2] = vec3{0.0F, 1.0F, depth};
    return face;
}

void testFiniteAndSafeMath() {
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float infinity = std::numeric_limits<float>::infinity();

    expect(isFinite(1.0F), "A finite scalar must be finite.");
    expect(!isFinite(nan), "NaN must not be finite.");
    expect(!isFinite(infinity), "Infinity must not be finite.");
    expect(isFinite(vec3{1.0F, 2.0F, 3.0F}), "A vector with finite components must be finite.");
    expect(!isFinite(vec3{1.0F, nan, 3.0F}), "Every vector component must be finite.");
    expect(!isFinite(vec2{infinity, 1.0F}), "A vector containing infinity must not be finite.");

    const vec3 normalized = safeNormalize(vec3{3.0F, 0.0F, 4.0F});
    expectNear(glm::length(normalized), 1.0F, 1.0e-6F, "safeNormalize must return a unit vector.");
    const vec3 fallback = safeNormalize(vec3{0.0F}, vec3{0.0F, 5.0F, 0.0F});
    expectNear(fallback.y, 1.0F, 1.0e-6F, "safeNormalize must normalize its fallback.");
    const vec3 invalidNormalized = safeNormalize(vec3{nan});
    expect(isFinite(invalidNormalized) && glm::length(invalidNormalized) == 0.0F,
           "safeNormalize must return zero when both the value and fallback are invalid.");

    expectNear(safeSqrt(-1.0F), 0.0F, 0.0F, "safeSqrt must clamp negative input.");
    expectNear(safeSqrt(nan), 0.0F, 0.0F, "safeSqrt must reject NaN input.");
    expectNear(safeSqrt(9.0F), 3.0F, 1.0e-6F, "safeSqrt must preserve valid square roots.");

    const vec3 powered = positivePow(vec3{-2.0F, 2.0F, 3.0F}, 2.0F);
    expectNear(powered.x, 0.0F, 0.0F, "positivePow must clamp negative bases.");
    expectNear(powered.y, 4.0F, 1.0e-6F, "positivePow must exponentiate positive bases.");
    expectNear(powered.z, 9.0F, 1.0e-6F, "positivePow must exponentiate every component.");
    expect(glm::length(positivePow(vec3{2.0F}, -1.0F)) == 0.0F,
           "positivePow must reject negative exponents without producing infinity.");
    expect(glm::length(positivePow(vec3{2.0F}, nan)) == 0.0F,
           "positivePow must reject non-finite exponents.");
}

void testGeometryIntersections() {
    const Box box{vec3{-1.0F, -1.0F, 2.0F}, vec3{1.0F, 1.0F, 4.0F}};
    const Ray ray{vec3{0.0F}, vec3{0.0F, 0.0F, 1.0F}};
    float nearDistance = 0.0F;
    float farDistance = std::numeric_limits<float>::infinity();
    expect(intersect(ray, box, nearDistance, farDistance), "A forward ray must intersect the box.");
    expectNear(nearDistance, 2.0F, 1.0e-6F, "The box entry distance must be correct.");
    expectNear(farDistance, 4.0F, 1.0e-6F, "The box exit distance must be correct.");

    nearDistance = 0.0F;
    farDistance = std::numeric_limits<float>::infinity();
    const Ray parallelMiss{vec3{2.0F, 0.0F, 0.0F}, vec3{0.0F, 0.0F, 1.0F}};
    expect(!intersect(parallelMiss, box, nearDistance, farDistance),
           "A parallel ray outside a slab must miss the box.");

    nearDistance = 0.0F;
    farDistance = std::numeric_limits<float>::infinity();
    expect(!intersect(ray, Box{}, nearDistance, farDistance), "An empty box must not intersect a ray.");

    const Face face = triangleAt(2.0F);
    expectNear(
        intersectTriangle(ray, face.vertices[0], face.vertices[1], face.vertices[2]),
        2.0F,
        1.0e-6F,
        "The ray-triangle distance must be correct.");
    expect(std::isinf(intersectTriangle(
               ray, vec3{0.0F}, vec3{1.0F, 0.0F, 0.0F}, vec3{2.0F, 0.0F, 0.0F})),
           "A degenerate triangle must not intersect.");

    const vec3 coordinates = barycentric(
        vec3{0.0F, 0.0F, 0.0F},
        vec3{1.0F, 0.0F, 0.0F},
        vec3{0.0F, 1.0F, 0.0F},
        vec3{0.25F, 0.25F, 0.0F});
    expectNear(coordinates.x, 0.5F, 1.0e-6F, "The first barycentric coordinate must be correct.");
    expectNear(coordinates.y, 0.25F, 1.0e-6F, "The second barycentric coordinate must be correct.");
    expectNear(coordinates.z, 0.25F, 1.0e-6F, "The third barycentric coordinate must be correct.");

    const vec3 degenerate = barycentric(
        vec3{0.0F}, vec3{1.0F, 0.0F, 0.0F}, vec3{2.0F, 0.0F, 0.0F}, vec3{0.5F});
    expect(isFinite(degenerate) && glm::length(degenerate) == 0.0F,
           "Degenerate barycentric coordinates must return a finite sentinel.");
}

void testTangentSpaces() {
    const vec3 normal{0.0F, 0.0F, 1.0F};
    vec3 tangent;
    vec3 bitangent;
    makeTangentSpace(normal, tangent, bitangent);
    expectUnitAndOrthogonal(normal, tangent, bitangent, "The generic tangent basis must be orthonormal.");

    makeTangentSpace(normal, vec3{0.0F, 0.0F, -1.0F}, tangent, bitangent);
    expectUnitAndOrthogonal(
        normal, tangent, bitangent, "A parallel incoming direction must use a valid fallback basis.");

    makeTangentSpace(vec3{0.0F}, tangent, bitangent);
    expectUnitAndOrthogonal(
        vec3{0.0F, 0.0F, 1.0F},
        tangent,
        bitangent,
        "A zero normal must use a valid fallback basis.");
}

void testRandomGeneratorAndDistribution() {
    Generator first{12345U};
    Generator second{12345U};
    for (int iteration = 0; iteration < 250'000; ++iteration) {
        const float value = first();
        expect(value > 0.0F && value < 1.0F, "Generator output must be strictly inside the unit interval.");
        if (iteration < 32) {
            expect(value == second(), "Generators with the same seed must be deterministic.");
        }
    }

    RandomDistribution distribution;
    Generator generator{7U};
    distribution.initialize({});
    expect(distribution.empty(), "An empty weight vector must create an empty distribution.");
    expect(distribution.sample(generator) == -1, "Sampling an empty distribution must fail safely.");
    expect(distribution.pdf(0) == 0.0F, "An empty distribution must report zero probability.");

    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float infinity = std::numeric_limits<float>::infinity();
    distribution.initialize({-1.0F, nan, 0.0F, infinity});
    expect(distribution.empty(), "A distribution with no positive finite weights must be empty.");
    expect(distribution.sample(generator) == -1, "Invalid weights must not produce an index.");

    distribution.initialize({-1.0F, nan, 0.0F, 1.0F, 3.0F, infinity});
    expect(!distribution.empty(), "Positive finite weights must remain sampleable.");
    expect(distribution.size() == 6, "Distribution indices must remain aligned with input weights.");
    expectNear(distribution.pdf(3), 0.25F, 1.0e-6F, "The first valid weight must have the correct PDF.");
    expectNear(distribution.pdf(4), 0.75F, 1.0e-6F, "The second valid weight must have the correct PDF.");
    expect(distribution.pdf(0) == 0.0F && distribution.pdf(1) == 0.0F &&
               distribution.pdf(2) == 0.0F && distribution.pdf(5) == 0.0F,
           "Invalid and non-positive weights must have zero probability.");
    expect(distribution.pdf(99) == 0.0F, "An out-of-range PDF query must fail safely.");

    for (int iteration = 0; iteration < 10'000; ++iteration) {
        const int selected = distribution.sample(generator);
        expect(selected == 3 || selected == 4, "Only positive finite weights may be sampled.");
    }
}

void testImageData() {
    ImageData image;
    expect(image.empty(), "A default image must be empty.");
    expect(!image.hasTransparency(), "An empty image must not report transparency.");
    expect(image.width() == 0 && image.height() == 0 && image.channels() == 0,
           "A default image must report zero dimensions and channels.");

    expectInvalidArgument(
        [&] { image.setBaseLevel(0, 1, 4, {}); },
        "Zero image dimensions must be rejected.");
    expectInvalidArgument(
        [&] { image.setBaseLevel(1, 1, 5, std::vector<std::uint8_t>(5)); },
        "Channel counts above four must be rejected.");
    expectInvalidArgument(
        [&] { image.setBaseLevel(2, 2, 4, std::vector<std::uint8_t>(15)); },
        "A mismatched image buffer must be rejected.");
    expectInvalidArgument(
        [&] {
            image.setBaseLevel(
                std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), 4, {});
        },
        "Oversized image dimensions must be rejected safely.");

    image.setBaseLevel(1, 1, 1, {64U});
    expect(image.width() == 1 && image.height() == 1 && image.channels() == 1,
           "Image metadata accessors must report the base level layout.");
    const vec4 gray = image.sampleRgba(0.25F, 0.75F, 0.0F);
    expectNear(gray.r, 64.0F / 255.0F, 1.0e-6F, "Grayscale sampling must populate red.");
    expectNear(gray.g, gray.r, 1.0e-6F, "Grayscale sampling must replicate green.");
    expectNear(gray.b, gray.r, 1.0e-6F, "Grayscale sampling must replicate blue.");
    expectNear(gray.a, 1.0F, 1.0e-6F, "Grayscale sampling must synthesize opaque alpha.");

    image.setBaseLevel(1, 1, 2, {128U, 17U});
    const vec4 grayAlpha = image.sampleRgba(0.0F, 0.0F, 0.0F);
    expectNear(grayAlpha.r, 128.0F / 255.0F, 1.0e-6F, "Gray-alpha sampling must preserve color.");
    expectNear(grayAlpha.a, 17.0F / 255.0F, 1.0e-6F, "Gray-alpha sampling must preserve alpha.");
    expect(image.hasTransparency(), "A sub-opaque alpha channel must report transparency.");

    image.setBaseLevel(1, 1, 3, {10U, 20U, 30U});
    const vec4 rgb = image.sampleRgba(0.0F, 0.0F, 0.0F);
    expectNear(rgb.r, 10.0F / 255.0F, 1.0e-6F, "RGB sampling must preserve red.");
    expectNear(rgb.g, 20.0F / 255.0F, 1.0e-6F, "RGB sampling must preserve green.");
    expectNear(rgb.b, 30.0F / 255.0F, 1.0e-6F, "RGB sampling must preserve blue.");
    expectNear(rgb.a, 1.0F, 1.0e-6F, "RGB sampling must synthesize opaque alpha.");
    expect(!image.hasTransparency(), "RGB data must not report transparency.");
    const float maximumCoordinate = std::numeric_limits<float>::max();
    const vec3 hugeCoordinateSample =
        image.sampleRgb(maximumCoordinate, -maximumCoordinate, maximumCoordinate);
    expectNear(hugeCoordinateSample.r, rgb.r, 1.0e-6F,
               "Huge finite texture coordinates must wrap before integer conversion.");
    expectNear(hugeCoordinateSample.g, rgb.g, 1.0e-6F,
               "Huge finite texture coordinates must preserve green.");
    expectNear(hugeCoordinateSample.b, rgb.b, 1.0e-6F,
               "Huge finite texture coordinates must preserve blue.");

    image.setBaseLevel(4, 4, 4, std::vector<std::uint8_t>(4U * 4U * 4U, 255U));
    image.generateMipmaps();
    expect(image.mipLevelCount() == 3, "A 4x4 image must produce 4x4, 2x2, and 1x1 mip levels.");
    const vec3 trilinear = image.sampleRgb(0.5F, 0.5F, 1.5F);
    expectNear(trilinear.r, 1.0F, 1.0e-6F, "Trilinear sampling must preserve constant colors.");
    expectNear(trilinear.g, 1.0F, 1.0e-6F, "Trilinear sampling must preserve constant colors.");
    expectNear(trilinear.b, 1.0F, 1.0e-6F, "Trilinear sampling must preserve constant colors.");

    image.setBaseLevel(3, 5, 1, std::vector<std::uint8_t>(3U * 5U, 64U));
    image.generateMipmaps();
    expect(image.mipLevelCount() == 4,
           "An odd non-square image must reduce through 3x5, 2x3, 1x2, and 1x1 levels.");

    image.setBaseLevel(3, 1, 1, {0U, 0U, 255U});
    image.generateMipmaps();
    const vec3 wrappedOddEdge = image.sampleRgb(0.5F, 0.5F, 1.0F);
    expectNear(wrappedOddEdge.r, 127.0F / 255.0F, 1.0e-6F,
               "Odd mip edges must wrap consistently with repeated texture sampling.");

    image.setBaseLevel(256, 1, 1, std::vector<std::uint8_t>(256U, 255U));
    image.generateMipmaps();
    expect(image.mipLevelCount() == 9,
           "Large textures must generate a complete mip chain instead of stopping after eight levels.");
}

void testMaterialConstants() {
    Material material;
    expect(!material.hasTexture(TextureSlot::count),
           "The texture-slot sentinel must not access material storage.");
    expectInvalidArgument(
        [&] { material.loadTexture(TextureSlot::count, {}); },
        "Loading the texture-slot sentinel must fail safely.");
    material.diffuseFactor = vec3{0.25F, 0.5F, 0.75F};
    material.emissiveFactor = vec3{2.0F, 1.0F, 0.5F};
    material.roughness = 0.35F;
    material.metallic = 0.65F;

    const vec4 diffuse = material.diffuseColor(0.0F, 0.0F, 0.0F);
    expectNear(diffuse.r, 0.25F, 1.0e-6F, "A textureless material must preserve constant diffuse red.");
    expectNear(diffuse.g, 0.5F, 1.0e-6F, "A textureless material must preserve constant diffuse green.");
    expectNear(diffuse.b, 0.75F, 1.0e-6F, "A textureless material must preserve constant diffuse blue.");
    expectNear(diffuse.a, 1.0F, 1.0e-6F, "A textureless material must remain opaque in its color value.");
    expect(material.isEmissive(), "A positive constant emission must identify an emissive material.");
    const vec3 emission = material.emissiveColor(0.0F, 0.0F, 0.0F);
    expectNear(emission.r, 2.0F, 1.0e-6F, "Constant emission must be returned without a texture.");

    float roughness = 0.0F;
    float metallic = 0.0F;
    material.surfaceParameters(0.0F, 0.0F, roughness, metallic);
    expectNear(roughness, 0.35F, 1.0e-6F, "Constant roughness must survive without a texture.");
    expectNear(metallic, 0.65F, 1.0e-6F, "Constant metallic must survive without a texture.");
}

void testPackedSceneValidation() {
    PackedSceneData scene;
    expectInvalidArgument(
        [&] { scene.validate(); }, "An empty packed scene must be rejected.");

    scene.vertices.push_back(PackedVertex{
        {0.0F, 0.0F, 1.0F, 0.0F},
        {0.0F, 1.0F, 0.25F, 0.75F},
    });
    scene.triangleIndices = {0U, 0U, 0U};
    scene.triangleMaterialIds = {0U};
    scene.materials.emplace_back();
    try {
        scene.validate();
    } catch (...) {
        expect(false, "A finite in-range packed scene must validate.");
    }

    scene.triangleIndices[2] = 1U;
    expectInvalidArgument(
        [&] { scene.validate(); }, "An out-of-range packed vertex index must be rejected.");
    scene.triangleIndices[2] = 0U;
    scene.vertices[0].normalYZAndUv[2] = std::numeric_limits<float>::quiet_NaN();
    expectInvalidArgument(
        [&] { scene.validate(); }, "A non-finite packed vertex attribute must be rejected.");
    scene.vertices[0].normalYZAndUv[2] = 0.25F;

    scene.materials[0].flagsAndReserved[0] = kPackedMaterialCutout;
    try {
        scene.validate();
    } catch (...) {
        expect(false, "A textureless cutout material must remain valid.");
    }
    scene.materials[0].flagsAndReserved[0] = 1U << 31U;
    expectInvalidArgument(
        [&] { scene.validate(); }, "An unknown packed material flag must be rejected.");
    scene.materials[0].flagsAndReserved[0] = 0U;

    scene.materials[0].textureIds[0] = 0U;
    expectInvalidArgument(
        [&] { scene.validate(); },
        "A texture ID without its matching presence flag must be rejected.");
    scene.materials[0].textureIds[0] = kInvalidSceneId;

    scene.materials[0].flagsAndReserved[1] = 1U;
    expectInvalidArgument(
        [&] { scene.validate(); }, "A nonzero reserved material flag lane must be rejected.");
    scene.materials[0].flagsAndReserved[1] = 0U;
    scene.materials[0].metallicSpecularAndReserved[2] = 1.0F;
    expectInvalidArgument(
        [&] { scene.validate(); }, "A nonzero reserved material float lane must be rejected.");
}

void testRenderSettingsValidation() {
    RenderSettings settings;
    settings.validate();
    expect(settings.resolvedThreadCount() >= 1, "Automatic worker selection must return at least one thread.");
    settings.threadCount = 100;
    settings.height = 3;
    expect(settings.resolvedThreadCount() == 3, "Worker selection must not exceed the image height.");

    settings = RenderSettings{};
    settings.width = 0;
    expectInvalidArgument([&] { settings.validate(); }, "Zero render width must be rejected.");
    settings = RenderSettings{};
    settings.samplesPerPixel = 0;
    expectInvalidArgument([&] { settings.validate(); }, "Zero samples per pixel must be rejected.");
    settings = RenderSettings{};
    settings.samplesPerPixel = 1'000'001;
    expectInvalidArgument([&] { settings.validate(); }, "Excessive sample counts must be rejected.");
    settings = RenderSettings{};
    settings.threadCount = -1;
    expectInvalidArgument([&] { settings.validate(); }, "Negative worker counts must be rejected.");
    settings = RenderSettings{};
    settings.direction.x = std::numeric_limits<float>::quiet_NaN();
    expectInvalidArgument([&] { settings.validate(); }, "Non-finite camera vectors must be rejected.");
    settings = RenderSettings{};
    settings.right = settings.direction;
    expectInvalidArgument([&] { settings.validate(); }, "A degenerate camera basis must be rejected.");
    settings = RenderSettings{};
    settings.directLightProbability = 1.1F;
    expectInvalidArgument([&] { settings.validate(); }, "Invalid direct-light probabilities must be rejected.");
    settings = RenderSettings{};
    settings.focusDistance = -1.0F;
    expectInvalidArgument([&] { settings.validate(); }, "Negative focus distances must be rejected.");
    settings = RenderSettings{};
    settings.exposure = 1.1e6F;
    expectInvalidArgument([&] { settings.validate(); }, "Excessive exposure must be rejected.");
    settings = RenderSettings{};
    settings.outputPrefix.clear();
    expectInvalidArgument([&] { settings.validate(); }, "An empty output prefix must be rejected.");
}

void testRadianceFinalization() {
    RadianceData radiance{vec3{1.0F, 2.0F, 3.0F}, 20.0F};
    finalizeRadianceData(radiance, 2.0F);
    expectNear(radiance.radiance.x, 2.0F, 0.0F, "Exposure must scale radiance X exactly.");
    expectNear(radiance.radiance.y, 4.0F, 0.0F, "Exposure must scale radiance Y exactly.");
    expectNear(radiance.radiance.z, 6.0F, 0.0F, "Exposure must scale radiance Z exactly.");
    expectNear(
        radiance.varianceAccumulator,
        24.0F,
        0.0F,
        "Radiance finalization must convert the exposed second moment to variance.");

    RadianceData invalid{
        vec3{std::numeric_limits<float>::quiet_NaN(), 1.0F, 2.0F},
        -1.0F};
    finalizeRadianceData(invalid, 1.0F);
    expect(
        invalid.radiance.x == 0.0F && invalid.radiance.y == 1.0F &&
            invalid.radiance.z == 2.0F &&
            invalid.varianceAccumulator == 0.0F,
        "Radiance finalization must sanitize NaN channels and invalid second moments.");
}

void testDielectricBsdf() {
    HitInfo surface;
    surface.shapeNormal = vec3{0.0F, 0.0F, 1.0F};
    surface.surfaceNormal = surface.shapeNormal;
    surface.position = vec3{0.0F};
    surface.baseColor = vec3{1.0F};
    surface.opacity = 0.0F;
    surface.roughness = 0.2F;
    surface.eta = 1.0F / 1.5F;
    surface.entering = true;

    Bsdf normalIncidence{vec3{0.0F, 0.0F, 1.0F}, surface};
    vec3 refracted;
    float fresnel = 0.0F;
    normalIncidence.refractPrecisely(refracted, fresnel);
    expect(isFinite(refracted) && refracted.z < 0.0F,
           "Normal-incidence dielectric refraction must enter the opposite hemisphere.");
    expectNear(fresnel, 0.04F, 1.0e-5F, "Air-to-glass normal-incidence Fresnel must be four percent.");

    const vec3 grazingIncoming = safeNormalize(vec3{0.9999F, 0.0F, 0.01F});
    Bsdf grazing{grazingIncoming, surface};
    grazing.refractPrecisely(refracted, fresnel);
    expect(isFinite(refracted), "Air-to-glass grazing refraction must remain finite.");
    expect(fresnel > 0.9F, "Dielectric Fresnel must approach one at grazing incidence.");

    surface.eta = 1.5F;
    surface.entering = false;
    const vec3 internalIncoming = safeNormalize(vec3{0.8660254F, 0.0F, 0.5F});
    Bsdf totalInternalReflection{internalIncoming, surface};
    totalInternalReflection.refractPrecisely(refracted, fresnel);
    expect(!isFinite(refracted) && fresnel == 1.0F,
           "Glass-to-air incidence above the critical angle must totally reflect.");

    Bsdf roughTransmission{vec3{0.0F, 0.0F, 1.0F}, surface};
    const vec3 transmittedValue = roughTransmission.evaluate(vec3{0.0F, 0.0F, -1.0F});
    expect(isFinite(transmittedValue) && transmittedValue.x > 0.0F &&
               transmittedValue.y > 0.0F && transmittedValue.z > 0.0F,
           "The rough dielectric BTDF must be finite and positive in a valid configuration.");

    surface.opacity = 1.0F;
    surface.entering = true;
    surface.baseColor = vec3{0.8F, 0.4F, 0.2F};
    Bsdf opaque{vec3{0.0F, 0.0F, 1.0F}, surface};
    const vec3 reflectedValue = opaque.evaluate(vec3{0.0F, 0.0F, 1.0F});
    expect(isFinite(reflectedValue) && reflectedValue.x >= 0.0F &&
               reflectedValue.y >= 0.0F && reflectedValue.z >= 0.0F,
           "The opaque reflection BSDF must remain finite and non-negative.");
}

void testBvh() {
    Bvh bvh;
    std::vector<Face> faces;
    bvh.build(faces);
    expect(bvh.empty(), "A BVH built from no faces must remain empty.");
    expect(bvh.nodeCount() == 0, "An empty BVH must not allocate nodes.");

    const Ray ray{vec3{0.0F}, vec3{0.0F, 0.0F, 1.0F}};
    HitRecord emptyHit;
    bvh.intersect(ray, emptyHit);
    expect(emptyHit.face == nullptr && std::isinf(emptyHit.tMaximum) &&
               emptyHit.primitiveIndex == std::numeric_limits<std::size_t>::max(),
           "Intersecting an empty BVH must preserve the miss record.");

    faces.push_back(triangleAt(2.0F));
    bvh.build(faces);
    expect(!bvh.empty() && bvh.nodeCount() == 1, "A one-face BVH must contain one leaf node.");
    HitRecord hit;
    bvh.intersect(ray, hit);
    expect(hit.face != nullptr, "The BVH must report a triangle hit.");
    expect(hit.primitiveIndex == 0U, "A BVH hit must report its face-array primitive index.");
    expectNear(hit.tMaximum, 2.0F, 1.0e-6F, "The BVH must report the nearest distance.");

    faces.clear();
    for (int depth = 25; depth >= 1; --depth) {
        faces.push_back(triangleAt(static_cast<float>(depth)));
    }
    bvh.build(faces);
    expect(bvh.nodeCount() > 1, "A larger BVH must contain internal nodes.");
    HitRecord nearestHit;
    bvh.intersect(ray, nearestHit);
    expect(nearestHit.face != nullptr, "A larger BVH must report a hit.");
    expect(nearestHit.primitiveIndex < faces.size() &&
               nearestHit.face == &faces[nearestHit.primitiveIndex],
           "A BVH hit primitive index must identify the returned face.");
    expectNear(nearestHit.tMaximum, 1.0F, 1.0e-6F, "BVH traversal must return the nearest triangle.");

    HitRecord clippedHit{1.5F, 3.5F};
    bvh.intersect(ray, clippedHit);
    expect(clippedHit.face != nullptr, "BVH traversal must honor a clipped query interval.");
    expectNear(
        clippedHit.tMaximum,
        2.0F,
        1.0e-6F,
        "BVH traversal must return the nearest triangle inside the query interval.");

    HitRecord boundedMiss{kRayEpsilon, 0.5F};
    bvh.intersect(ray, boundedMiss);
    expect(boundedMiss.face == nullptr && boundedMiss.tMaximum == 0.5F &&
               boundedMiss.primitiveIndex == std::numeric_limits<std::size_t>::max(),
           "BVH traversal must preserve a hit bound that ends before the root.");

    const Ray reverseRay{vec3{0.0F, 0.0F, 30.0F}, vec3{0.0F, 0.0F, -1.0F}};
    HitRecord reverseHit;
    bvh.intersect(reverseRay, reverseHit);
    expect(reverseHit.face != nullptr, "BVH traversal must support negative ray directions.");
    expectNear(
        reverseHit.tMaximum,
        5.0F,
        1.0e-6F,
        "BVH traversal must visit the geometrically nearer child first for a reverse ray.");

    faces.clear();
    bvh.build(faces);
    expect(bvh.empty() && bvh.nodeCount() == 0, "Rebuilding with no faces must clear the BVH.");
}

}  // namespace

int main() {
    testFiniteAndSafeMath();
    testGeometryIntersections();
    testTangentSpaces();
    testRandomGeneratorAndDistribution();
    testImageData();
    testMaterialConstants();
    testPackedSceneValidation();
    testRenderSettingsValidation();
    testRadianceFinalization();
    testDielectricBsdf();
    testBvh();

    if (failureCount != 0) {
        std::cerr << failureCount << " core test assertion(s) failed.\n";
        return 1;
    }

    std::cout << "All core tests passed.\n";
    return 0;
}
