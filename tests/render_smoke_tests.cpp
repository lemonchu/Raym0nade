#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "raym0nade/image.hpp"
#include "raym0nade/model.hpp"
#include "raym0nade/render.hpp"

#ifndef RAYM0NADE_TEST_SOURCE_DIR
#error "RAYM0NADE_TEST_SOURCE_DIR must be defined by CMake."
#endif

#ifndef RAYM0NADE_TEST_BINARY_DIR
#error "RAYM0NADE_TEST_BINARY_DIR must be defined by CMake."
#endif

namespace {

using namespace raym0nade;

constexpr std::array<const char*, 20> kOutputTags{
    "DiffuseColor",
    "DiffuseColor_FXAA",
    "ShapeNormal",
    "SurfaceNormal",
    "Direct_Diffuse",
    "Direct_Specular",
    "Indirect_Diffuse",
    "Indirect_Specular",
    "Raw",
    "Raw_Bloom",
    "Raw_FXAA",
    "Raw_Bloom_FXAA",
    "Direct_Diffuse_Filter",
    "Direct_Specular_Filter",
    "Indirect_Diffuse_Filter",
    "Indirect_Specular_Filter",
    "Filter",
    "Filter_Bloom",
    "Filter_FXAA",
    "Filter_Bloom_FXAA",
};

std::vector<std::uint8_t> readBytes(const std::filesystem::path& filename) {
    std::ifstream stream{filename, std::ios::binary};
    if (!stream) {
        throw std::runtime_error("Could not open test output: " + filename.string());
    }
    return std::vector<std::uint8_t>{
        std::istreambuf_iterator<char>{stream}, std::istreambuf_iterator<char>{}};
}

std::filesystem::path outputFile(
    const std::filesystem::path& prefix, const std::string& tag) {
    std::filesystem::path filename = prefix;
    filename += "(" + tag + ").png";
    return filename;
}

}  // namespace

int main() {
    try {
        const std::filesystem::path sourceDirectory{RAYM0NADE_TEST_SOURCE_DIR};
        const std::filesystem::path outputDirectory =
            std::filesystem::path{RAYM0NADE_TEST_BINARY_DIR} / "render-smoke-output";
        std::filesystem::create_directories(outputDirectory);

        Model invalidVertexModel{
            sourceDirectory / "tests" / "data", "non_finite_vertices.obj", "null"};
        if (invalidVertexModel.faceCount() != 1) {
            throw std::runtime_error(
                "Faces with non-finite vertex coordinates must be skipped during import.");
        }
        const HitRecord validTriangleHit =
            invalidVertexModel.intersect(Ray{vec3{0.0F}, vec3{0.0F, 0.0F, 1.0F}});
        if (validTriangleHit.face == nullptr) {
            throw std::runtime_error(
                "A rejected non-finite face must not hide valid imported geometry.");
        }

        Model model{sourceDirectory / "tests" / "data", "triangle.obj", "null"};
        if (model.faceCount() != 2 || model.lights().empty()) {
            throw std::runtime_error(
                "The smoke-test model must contain one diffuse face and one area light.");
        }
        const PackedSceneData packedScene = model.packScene();
        if (packedScene.triangleCount() != model.faceCount() ||
            packedScene.triangleMaterialIds.size() != model.faceCount() ||
            packedScene.vertices.empty() || packedScene.materials.empty() ||
            packedScene.areaLights.size() != 1U ||
            packedScene.areaLightTriangles.size() != 1U) {
            throw std::runtime_error(
                "The imported smoke scene did not produce complete packed geometry and lighting.");
        }
        packedScene.validate();
        const PackedAreaLight& packedLight = packedScene.areaLights.front();
        const PackedAreaLightTriangle& packedLightTriangle =
            packedScene.areaLightTriangles.front();
        if (packedLight.centerAndPower[3] <= 0.0F ||
            packedLight.triangleRangeAndReserved[0] != 0U ||
            packedLight.triangleRangeAndReserved[1] != 1U ||
            packedLightTriangle.areaProbabilityCdfAndReserved[0] <= 0.0F ||
            packedLightTriangle.areaProbabilityCdfAndReserved[1] != 1.0F ||
            packedLightTriangle.areaProbabilityCdfAndReserved[2] != 1.0F) {
            throw std::runtime_error(
                "The imported smoke-scene area light has invalid packed sampling data.");
        }
        const HitRecord centerHit = model.intersect(Ray{vec3{0.0F}, vec3{0.0F, 0.0F, 1.0F}});
        if (centerHit.face == nullptr) {
            throw std::runtime_error("The smoke-test ray did not hit its triangle.");
        }
        if (centerHit.primitiveIndex >= packedScene.triangleCount()) {
            throw std::runtime_error(
                "The CPU hit did not map to the packed-scene primitive range.");
        }

        RenderSettings settings;
        settings.width = 8;
        settings.height = 8;
        settings.samplesPerPixel = 1;
        settings.directLightProbability = 0.5F;
        settings.pixelScale = 0.125F;
        settings.seed = 0xC0FFEEU;
        settings.threadCount = 1;

        const std::filesystem::path filmOnlyDirectory = outputDirectory / "film-only";
        std::filesystem::remove_all(filmOnlyDirectory);
        settings.outputPrefix = filmOnlyDirectory / "must-not-be-written";
        FilmRenderResult filmOnlyResult = renderToFilm(model, settings);
        if (std::filesystem::exists(filmOnlyDirectory)) {
            throw std::runtime_error(
                "Rendering to a Film must not create output directories or files.");
        }
        if (filmOnlyResult.film.width() != settings.width ||
            filmOnlyResult.film.height() != settings.height ||
            filmOnlyResult.stats.renderSeconds < 0.0 ||
            filmOnlyResult.stats.totalSeconds != filmOnlyResult.stats.renderSeconds ||
            filmOnlyResult.stats.directLightSamples == 0) {
            throw std::runtime_error(
                "The in-memory Film render did not return complete pixels and statistics.");
        }

        const std::filesystem::path singleThreadPrefix = outputDirectory / "single";
        settings.outputPrefix = singleThreadPrefix;
        const RenderStats singleThreadStats = renderToFiles(model, settings);

        settings.threadCount = 4;
        const std::filesystem::path multiThreadPrefix = outputDirectory / "multi";
        settings.outputPrefix = multiThreadPrefix;
        const RenderStats multiThreadStats = renderToFiles(model, settings);
        if (singleThreadStats.renderSeconds < 0.0 || multiThreadStats.renderSeconds < 0.0) {
            throw std::runtime_error("Render timing must never be negative.");
        }
        if (singleThreadStats.directLightSamples == 0 ||
            singleThreadStats.directLightSamples != filmOnlyResult.stats.directLightSamples ||
            singleThreadStats.directLightSamples != multiThreadStats.directLightSamples) {
            throw std::runtime_error(
                "The smoke test must exercise deterministic random direct-light sampling.");
        }

        for (const char* tag : kOutputTags) {
            const std::vector<std::uint8_t> single = readBytes(outputFile(singleThreadPrefix, tag));
            const std::vector<std::uint8_t> multi = readBytes(outputFile(multiThreadPrefix, tag));
            if (single.empty() || single != multi) {
                throw std::runtime_error(
                    std::string{"Worker-count determinism failed for output: "} + tag);
            }
        }

        Film diffuse{outputFile(singleThreadPrefix, "DiffuseColor")};
        bool foundMaterialColor = false;
        for (const vec3& pixel : diffuse.pixels) {
            foundMaterialColor = foundMaterialColor ||
                                 (isFinite(pixel) && pixel.r > pixel.g && pixel.g > pixel.b &&
                                  pixel.r > 0.1F);
        }
        if (!foundMaterialColor) {
            throw std::runtime_error(
                "The imported constant diffuse material was not visible in the diagnostic pass.");
        }

        Film directDiffuse{outputFile(singleThreadPrefix, "Direct_Diffuse")};
        bool foundDirectLighting = false;
        for (const vec3& pixel : directDiffuse.pixels) {
            foundDirectLighting = foundDirectLighting ||
                                  (isFinite(pixel) && glm::dot(pixel, pixel) > 1.0e-6F);
        }
        if (!foundDirectLighting) {
            throw std::runtime_error(
                "The area light did not produce a visible direct-light sample.");
        }

        const std::filesystem::path linearImagePath =
            outputDirectory / "linear-image-display-scale.png";
        LinearImage linearImage{
            ImageExtent{2U, 1U},
            {
                vec3{std::numeric_limits<float>::max()},
                vec3{0.25F, 0.125F, 0.0625F},
            },
        };
        saveLinearImagePng(linearImage, linearImagePath, 2.0F);
        const Film scaledImage{linearImagePath};
        if (scaledImage.width() != 2 || scaledImage.height() != 1 ||
            scaledImage.pixels[0].x < 0.9F || scaledImage.pixels[0].y < 0.9F ||
            scaledImage.pixels[0].z < 0.9F || scaledImage.pixels[1].x <= 0.0F) {
            throw std::runtime_error(
                "Linear-image display scaling must keep overflowing highlights bright.");
        }
        bool rejectedInvalidScale = false;
        try {
            saveLinearImagePng(linearImage, linearImagePath, -1.0F);
        } catch (const std::invalid_argument&) {
            rejectedInvalidScale = true;
        }
        if (!rejectedInvalidScale) {
            throw std::runtime_error(
                "Linear-image PNG export must reject a negative display scale.");
        }

        std::cout << "Render smoke test passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
