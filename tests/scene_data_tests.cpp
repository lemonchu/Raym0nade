#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "raym0nade/image.hpp"
#include "raym0nade/model.hpp"
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

template <typename Function>
void expectInvalidArgument(Function function, std::string_view message) {
    try {
        function();
        expect(false, message);
    } catch (const std::invalid_argument&) {
    } catch (...) {
        expect(false, message);
    }
}

class TemporaryDirectory {
public:
    TemporaryDirectory() {
        const std::filesystem::path base = std::filesystem::temp_directory_path();
        const auto stamp =
            std::chrono::steady_clock::now().time_since_epoch().count();
        for (unsigned int attempt = 0U; attempt < 64U; ++attempt) {
            std::error_code error;
            const std::filesystem::path candidate =
                base / ("raym0nade-scene-data-" + std::to_string(stamp) + "-" +
                        std::to_string(attempt));
            if (std::filesystem::create_directory(candidate, error)) {
                path_ = candidate;
                return;
            }
        }
        throw std::runtime_error("Could not create the scene-data test directory.");
    }

    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

    ~TemporaryDirectory() {
        std::error_code error;
        if (!path_.empty()) {
            std::filesystem::remove_all(path_, error);
        }
    }

    [[nodiscard]] const std::filesystem::path& path() const noexcept {
        return path_;
    }

private:
    std::filesystem::path path_;
};

void writeTextFile(const std::filesystem::path& path, std::string_view contents) {
    std::ofstream output{path, std::ios::binary};
    if (!output) {
        throw std::runtime_error("Could not create test fixture: " + path.string());
    }
    output << contents;
    if (!output) {
        throw std::runtime_error("Could not write test fixture: " + path.string());
    }
}

void writeTinyHdr(const std::filesystem::path& path) {
    std::ofstream output{path, std::ios::binary};
    if (!output) {
        throw std::runtime_error("Could not create the HDR test fixture.");
    }
    constexpr std::string_view header{
        "#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y 1 +X 2\n"};
    constexpr std::array<std::uint8_t, 8> pixels{
        128U, 0U, 0U, 129U,
        0U, 128U, 0U, 129U,
    };
    output.write(header.data(), static_cast<std::streamsize>(header.size()));
    output.write(
        reinterpret_cast<const char*>(pixels.data()),
        static_cast<std::streamsize>(pixels.size()));
    if (!output) {
        throw std::runtime_error("Could not write the HDR test fixture.");
    }
}

PackedSceneData makeGeometryScene() {
    PackedSceneData scene;
    scene.vertices = {
        PackedVertex{
            {0.0F, 0.0F, 1.0F, 0.0F},
            {0.0F, -1.0F, 0.0F, 0.0F},
        },
        PackedVertex{
            {1.0F, 0.0F, 1.0F, 0.0F},
            {0.0F, -1.0F, 1.0F, 0.0F},
        },
        PackedVertex{
            {0.0F, 1.0F, 1.0F, 0.0F},
            {0.0F, -1.0F, 0.0F, 1.0F},
        },
    };
    scene.triangleIndices = {0U, 1U, 2U};
    scene.triangleMaterialIds = {0U};
    scene.materials.emplace_back();
    return scene;
}

void addOnePixelDiffuseTexture(PackedSceneData& scene) {
    scene.textures.push_back(PackedTexture{0U, 1U, 1U, 1U});
    scene.textureMipLevels.push_back(PackedTextureMip{0U, 1U, 1U, 1U});
    scene.textureTexelsRgba8.push_back(0xff332211U);
    scene.materials[0].flagsAndReserved[0] |= kPackedMaterialHasDiffuseTexture;
    scene.materials[0].textureIds[0] = 0U;
}

PackedSceneData makeOddMipScene() {
    PackedSceneData scene = makeGeometryScene();
    scene.textures.push_back(PackedTexture{0U, 4U, 3U, 5U});
    scene.textureMipLevels = {
        PackedTextureMip{0U, 15U, 3U, 5U},
        PackedTextureMip{15U, 6U, 2U, 3U},
        PackedTextureMip{21U, 2U, 1U, 2U},
        PackedTextureMip{23U, 1U, 1U, 1U},
    };
    scene.textureTexelsRgba8.assign(24U, 0xffffffffU);
    scene.materials[0].flagsAndReserved[0] |= kPackedMaterialHasDiffuseTexture;
    scene.materials[0].textureIds[0] = 0U;
    return scene;
}

void addOneAreaLight(PackedSceneData& scene) {
    scene.materials[0].emissionAndIor[0] = 2.0F;
    scene.areaLights.push_back(PackedAreaLight{
        {0.0F, 0.0F, 1.0F, 2.0F},
        {0U, 1U, 0U, 0U},
    });
    scene.areaLightTriangles.push_back(PackedAreaLightTriangle{
        {0U, 1U, 2U, 0U},
        {0.5F, 1.0F, 1.0F, 0.0F},
    });
}

void addSampleableEnvironment(PackedSceneData& scene) {
    scene.environment =
        PackedEnvironment{2U, 2U, kPackedEnvironmentHasImportance, 0U};
    scene.environmentRows = {
        PackedEnvironmentRow{3.0F / 11.0F, 3.0F / 11.0F, kPi, 0.0F},
        PackedEnvironmentRow{8.0F / 11.0F, 1.0F, kPi, 0.0F},
    };
    scene.environmentTexels = {
        PackedEnvironmentTexel{{1.0F, 0.0F, 0.0F, 1.0F / 3.0F}},
        PackedEnvironmentTexel{{2.0F, 0.0F, 0.0F, 1.0F}},
        PackedEnvironmentTexel{{0.0F, 1.0F, 0.0F, 0.25F}},
        PackedEnvironmentTexel{{0.0F, 3.0F, 0.0F, 1.0F}},
    };
}

void testLightingValidation() {
    PackedSceneData lightScene = makeGeometryScene();
    addOneAreaLight(lightScene);
    try {
        lightScene.validate();
    } catch (...) {
        expect(false, "A complete packed area light must validate.");
    }

    PackedSceneData badLightRange = lightScene;
    badLightRange.areaLights[0].triangleRangeAndReserved[0] = 1U;
    expectInvalidArgument(
        [&] { badLightRange.validate(); },
        "A non-contiguous area-light range must be rejected.");

    PackedSceneData orphanLightTriangle = lightScene;
    orphanLightTriangle.areaLights.clear();
    expectInvalidArgument(
        [&] { orphanLightTriangle.validate(); },
        "An unreferenced area-light triangle must be rejected.");

    PackedSceneData badLightVertex = lightScene;
    badLightVertex.areaLightTriangles[0].vertexIdsAndMaterialId[0] = 1U;
    expectInvalidArgument(
        [&] { badLightVertex.validate(); },
        "An out-of-range area-light vertex must be rejected.");

    PackedSceneData nonEmissiveLight = lightScene;
    nonEmissiveLight.materials[0].emissionAndIor[0] = 0.0F;
    expectInvalidArgument(
        [&] { nonEmissiveLight.validate(); },
        "An area-light triangle with a non-emissive material must be rejected.");

    PackedSceneData badFaceDistribution = lightScene;
    badFaceDistribution.areaLightTriangles[0]
        .areaProbabilityCdfAndReserved[2] = 0.5F;
    expectInvalidArgument(
        [&] { badFaceDistribution.validate(); },
        "An area-light face CDF that does not end at one must be rejected.");

    PackedSceneData badLightArea = lightScene;
    badLightArea.areaLightTriangles[0].areaProbabilityCdfAndReserved[0] =
        0.25F;
    expectInvalidArgument(
        [&] { badLightArea.validate(); },
        "An area-light area that disagrees with its triangle must be rejected.");

    PackedSceneData environmentScene = makeGeometryScene();
    addSampleableEnvironment(environmentScene);
    try {
        environmentScene.validate();
    } catch (...) {
        expect(false, "A complete sampleable packed environment must validate.");
    }

    PackedSceneData blackEnvironment = makeGeometryScene();
    blackEnvironment.environment = PackedEnvironment{2U, 1U, 0U, 0U};
    blackEnvironment.environmentRows = {
        PackedEnvironmentRow{0.0F, 0.0F, 2.0F * kPi, 0.0F},
    };
    blackEnvironment.environmentTexels = {
        PackedEnvironmentTexel{{0.0F, 0.0F, 0.0F, 0.5F}},
        PackedEnvironmentTexel{{0.0F, 0.0F, 0.0F, 1.0F}},
    };
    try {
        blackEnvironment.validate();
    } catch (...) {
        expect(false, "A black lookup-only packed environment must validate.");
    }

    PackedSceneData zeroWeightRow = makeGeometryScene();
    zeroWeightRow.environment =
        PackedEnvironment{1U, 2U, kPackedEnvironmentHasImportance, 0U};
    zeroWeightRow.environmentRows = {
        PackedEnvironmentRow{1.0F, 1.0F, 2.0F * kPi, 0.0F},
        PackedEnvironmentRow{0.0F, 1.0F, 2.0F * kPi, 0.0F},
    };
    zeroWeightRow.environmentTexels = {
        PackedEnvironmentTexel{{1.0F, 0.0F, 0.0F, 1.0F}},
        PackedEnvironmentTexel{{0.0F, 0.0F, 0.0F, 1.0F}},
    };
    try {
        zeroWeightRow.validate();
    } catch (...) {
        expect(
            false,
            "A zero-weight row in a sampleable environment must retain a valid CDF plateau.");
    }

    PackedSceneData badEnvironmentExtent = environmentScene;
    badEnvironmentExtent.environment.width = 0U;
    expectInvalidArgument(
        [&] { badEnvironmentExtent.validate(); },
        "A partially zero environment extent must be rejected.");

    PackedSceneData badEnvironmentFlag = environmentScene;
    badEnvironmentFlag.environment.flags |= 0x80000000U;
    expectInvalidArgument(
        [&] { badEnvironmentFlag.validate(); },
        "An unknown packed-environment flag must be rejected.");

    PackedSceneData missingEnvironmentTexel = environmentScene;
    missingEnvironmentTexel.environmentTexels.pop_back();
    expectInvalidArgument(
        [&] { missingEnvironmentTexel.validate(); },
        "Environment storage that does not match its extent must be rejected.");

    PackedSceneData badRowDistribution = environmentScene;
    badRowDistribution.environmentRows[0].probability = 0.5F;
    expectInvalidArgument(
        [&] { badRowDistribution.validate(); },
        "Inconsistent environment row probability and CDF data must be rejected.");

    PackedSceneData badConditionalCdf = environmentScene;
    badConditionalCdf.environmentTexels[1].radianceAndConditionalCdf[3] = 0.9F;
    expectInvalidArgument(
        [&] { badConditionalCdf.validate(); },
        "An environment conditional CDF that does not end at one must be rejected.");

    PackedSceneData biasedConditionalCdf = environmentScene;
    biasedConditionalCdf.environmentTexels[0].radianceAndConditionalCdf[3] =
        0.4F;
    expectInvalidArgument(
        [&] { biasedConditionalCdf.validate(); },
        "An environment CDF that disagrees with its radiance must be rejected.");

    PackedSceneData badSolidAngle = environmentScene;
    badSolidAngle.environmentRows[0].solidAngle = 1.0F;
    expectInvalidArgument(
        [&] { badSolidAngle.validate(); },
        "An environment solid angle that disagrees with its extent must be rejected.");

    PackedSceneData missingImportance = environmentScene;
    missingImportance.environment.flags = 0U;
    for (PackedEnvironmentRow& row : missingImportance.environmentRows) {
        row.probability = 0.0F;
        row.cumulativeProbability = 0.0F;
    }
    expectInvalidArgument(
        [&] { missingImportance.validate(); },
        "A non-black environment without importance metadata must be rejected.");

    PackedSceneData spuriousImportance = blackEnvironment;
    spuriousImportance.environment.flags = kPackedEnvironmentHasImportance;
    spuriousImportance.environmentRows[0].probability = 1.0F;
    spuriousImportance.environmentRows[0].cumulativeProbability = 1.0F;
    expectInvalidArgument(
        [&] { spuriousImportance.validate(); },
        "A black environment with importance metadata must be rejected.");

    PackedSceneData negativeEnvironment = environmentScene;
    negativeEnvironment.environmentTexels[0].radianceAndConditionalCdf[0] = -1.0F;
    expectInvalidArgument(
        [&] { negativeEnvironment.validate(); },
        "Negative packed-environment radiance must be rejected.");
}

void testTextureValidation() {
    PackedSceneData untextured = makeGeometryScene();
    try {
        untextured.validate();
    } catch (...) {
        expect(false, "A scene with only textureless materials must remain valid.");
    }

    PackedSceneData onePixel = makeGeometryScene();
    addOnePixelDiffuseTexture(onePixel);
    try {
        onePixel.validate();
    } catch (...) {
        expect(false, "A complete one-pixel texture store must validate.");
    }

    PackedSceneData flagWithoutId = untextured;
    flagWithoutId.materials[0].flagsAndReserved[0] |=
        kPackedMaterialHasDiffuseTexture;
    expectInvalidArgument(
        [&] { flagWithoutId.validate(); },
        "A texture presence flag without an ID must be rejected.");

    PackedSceneData idWithoutFlag = untextured;
    idWithoutFlag.materials[0].textureIds[0] = 0U;
    expectInvalidArgument(
        [&] { idWithoutFlag.validate(); },
        "A texture ID without a presence flag must be rejected.");

    PackedSceneData outOfRangeId = onePixel;
    outOfRangeId.materials[0].textureIds[0] = 1U;
    expectInvalidArgument(
        [&] { outOfRangeId.validate(); },
        "An out-of-range texture ID must be rejected.");

    PackedSceneData zeroExtent = onePixel;
    zeroExtent.textures[0].width = 0U;
    expectInvalidArgument(
        [&] { zeroExtent.validate(); },
        "A zero texture extent must be rejected.");

    PackedSceneData badMipRange = onePixel;
    badMipRange.textures[0].firstMipLevel = 1U;
    expectInvalidArgument(
        [&] { badMipRange.validate(); },
        "An out-of-range texture mip range must be rejected.");

    PackedSceneData badTexelOffset = onePixel;
    badTexelOffset.textureMipLevels[0].texelOffset = 1U;
    expectInvalidArgument(
        [&] { badTexelOffset.validate(); },
        "A non-contiguous texture texel offset must be rejected.");

    PackedSceneData badTexelCount = onePixel;
    badTexelCount.textureMipLevels[0].texelCount = 2U;
    expectInvalidArgument(
        [&] { badTexelCount.validate(); },
        "A mip texel count that disagrees with its extent must be rejected.");

    PackedSceneData orphanMip = onePixel;
    orphanMip.textureMipLevels.push_back(PackedTextureMip{1U, 1U, 1U, 1U});
    expectInvalidArgument(
        [&] { orphanMip.validate(); },
        "An unreferenced mip descriptor must be rejected.");

    PackedSceneData orphanTexel = onePixel;
    orphanTexel.textureTexelsRgba8.push_back(0U);
    expectInvalidArgument(
        [&] { orphanTexel.validate(); },
        "An unreferenced packed texel must be rejected.");

    PackedSceneData repeatedTerminalMip = onePixel;
    repeatedTerminalMip.textures[0].mipLevelCount = 2U;
    repeatedTerminalMip.textureMipLevels.push_back(
        PackedTextureMip{1U, 1U, 1U, 1U});
    repeatedTerminalMip.textureTexelsRgba8.push_back(0xff332211U);
    expectInvalidArgument(
        [&] { repeatedTerminalMip.validate(); },
        "A repeated terminal 1x1 mip must be rejected.");

    PackedSceneData incomplete = makeGeometryScene();
    incomplete.textures.push_back(PackedTexture{0U, 1U, 3U, 5U});
    incomplete.textureMipLevels.push_back(PackedTextureMip{0U, 15U, 3U, 5U});
    incomplete.textureTexelsRgba8.assign(15U, 0U);
    incomplete.materials[0].flagsAndReserved[0] |=
        kPackedMaterialHasDiffuseTexture;
    incomplete.materials[0].textureIds[0] = 0U;
    expectInvalidArgument(
        [&] { incomplete.validate(); },
        "An incomplete ceil-halved mip chain must be rejected.");

    PackedSceneData wrongOddMip = makeOddMipScene();
    wrongOddMip.textureMipLevels[1].height = 2U;
    expectInvalidArgument(
        [&] { wrongOddMip.validate(); },
        "A malformed odd-sized mip chain must be rejected.");

    PackedSceneData overflowing = onePixel;
    overflowing.textures[0].width = std::numeric_limits<std::uint32_t>::max();
    overflowing.textures[0].height = std::numeric_limits<std::uint32_t>::max();
    overflowing.textureMipLevels[0].width =
        std::numeric_limits<std::uint32_t>::max();
    overflowing.textureMipLevels[0].height =
        std::numeric_limits<std::uint32_t>::max();
    expectInvalidArgument(
        [&] { overflowing.validate(); },
        "Overflowing texture dimensions must be rejected without allocation.");

    PackedSceneData oddMip = makeOddMipScene();
    try {
        oddMip.validate();
    } catch (...) {
        expect(false, "A complete 3x5 ceil-halved mip chain must validate.");
    }
}

void createImportedFixture(const std::filesystem::path& directory) {
    const std::filesystem::path textureDirectory = directory / "textures";
    std::filesystem::create_directories(textureDirectory);

    Film oddTexture{3, 5};
    std::fill(
        oddTexture.pixels.begin(), oddTexture.pixels.end(), vec3{1.0F, 0.0F, 0.0F});
    oddTexture.save(textureDirectory / "odd.png");

    Film uniqueTexture{2, 1};
    std::fill(
        uniqueTexture.pixels.begin(),
        uniqueTexture.pixels.end(),
        vec3{0.0F, 1.0F, 0.0F});
    uniqueTexture.save(textureDirectory / "unique.png");

    writeTextFile(
        directory / "scene.mtl",
        "newmtl SharedFirst\n"
        "Kd 1 1 1\n"
        "map_Kd textures/odd.png\n"
        "\n"
        "newmtl SharedSecond\n"
        "Kd 1 1 1\n"
        "map_Kd textures/./odd.png\n"
        "\n"
        "newmtl Unique\n"
        "Kd 1 1 1\n"
        "map_Kd textures/unique.png\n");

    writeTextFile(
        directory / "scene.obj",
        "mtllib scene.mtl\n"
        "o PackedTextures\n"
        "v -3 -1 3\n"
        "v -1 -1 3\n"
        "v -2 1 3\n"
        "v -1 -1 3\n"
        "v 1 -1 3\n"
        "v 0 1 3\n"
        "v 1 -1 3\n"
        "v 3 -1 3\n"
        "v 2 1 3\n"
        "vt 0 0\n"
        "vt 1 0\n"
        "vt 0.5 1\n"
        "vn 0 0 -1\n"
        "usemtl SharedFirst\n"
        "f 1/1/1 2/2/1 3/3/1\n"
        "usemtl SharedSecond\n"
        "f 4/1/1 5/2/1 6/3/1\n"
        "usemtl Unique\n"
        "f 7/1/1 8/2/1 9/3/1\n");
}

bool sameTexture(const PackedTexture& left, const PackedTexture& right) {
    return left.firstMipLevel == right.firstMipLevel &&
           left.mipLevelCount == right.mipLevelCount &&
           left.width == right.width && left.height == right.height;
}

bool sameMip(const PackedTextureMip& left, const PackedTextureMip& right) {
    return left.texelOffset == right.texelOffset &&
           left.texelCount == right.texelCount && left.width == right.width &&
           left.height == right.height;
}

void testImportedTexturePacking() {
    TemporaryDirectory directory;
    createImportedFixture(directory.path());

    Model model{directory.path(), "scene.obj", "null"};
    const PackedSceneData first = model.packScene();
    const PackedSceneData second = model.packScene();
    first.validate();
    second.validate();

    expect(first.textures.size() == 2U,
           "Normalized source paths must deduplicate shared textures.");
    expect(first.textureMipLevels.size() == 6U,
           "The two imported textures must retain all generated mip levels.");
    expect(first.textureTexelsRgba8.size() == 27U,
           "Packed mip ranges must occupy exactly 27 RGBA8 texels.");

    if (first.textures.size() == 2U && first.textureMipLevels.size() == 6U) {
        const PackedTexture& odd = first.textures[0];
        expect(
            odd.firstMipLevel == 0U && odd.mipLevelCount == 4U &&
                odd.width == 3U && odd.height == 5U,
            "The first stable texture ID must describe the 3x5 source.");
        const std::array<PackedTextureMip, 4> expectedOddMips{
            PackedTextureMip{0U, 15U, 3U, 5U},
            PackedTextureMip{15U, 6U, 2U, 3U},
            PackedTextureMip{21U, 2U, 1U, 2U},
            PackedTextureMip{23U, 1U, 1U, 1U},
        };
        for (std::size_t level = 0U; level < expectedOddMips.size(); ++level) {
            expect(
                sameMip(first.textureMipLevels[level], expectedOddMips[level]),
                "Odd-sized mip metadata must preserve ceil-halved offsets.");
        }

        const PackedTexture& unique = first.textures[1];
        expect(
            unique.firstMipLevel == 4U && unique.mipLevelCount == 2U &&
                unique.width == 2U && unique.height == 1U,
            "The second stable texture ID must begin after the odd mip chain.");
        expect(
            sameMip(
                first.textureMipLevels[4],
                PackedTextureMip{24U, 2U, 2U, 1U}) &&
                sameMip(
                    first.textureMipLevels[5],
                    PackedTextureMip{26U, 1U, 1U, 1U}),
            "The second texture must use contiguous global texel offsets.");
    }

    constexpr std::uint32_t kOpaqueRedRgba8 = 0xff0000ffU;
    constexpr std::uint32_t kOpaqueGreenRgba8 = 0xff00ff00U;
    if (first.textureTexelsRgba8.size() == 27U) {
        expect(
            std::all_of(
                first.textureTexelsRgba8.begin(),
                first.textureTexelsRgba8.begin() + 24,
                [](std::uint32_t texel) { return texel == kOpaqueRedRgba8; }),
            "Encoded red bytes and synthesized alpha must survive RGBA8 packing.");
        expect(
            std::all_of(
                first.textureTexelsRgba8.begin() + 24,
                first.textureTexelsRgba8.end(),
                [](std::uint32_t texel) { return texel == kOpaqueGreenRgba8; }),
            "Encoded green bytes and synthesized alpha must survive RGBA8 packing.");
    }

    std::vector<std::uint32_t> diffuseTextureIds;
    for (const PackedMaterial& material : first.materials) {
        if ((material.flagsAndReserved[0] &
             kPackedMaterialHasDiffuseTexture) != 0U) {
            diffuseTextureIds.push_back(material.textureIds[0]);
        }
    }
    expect(
        diffuseTextureIds == std::vector<std::uint32_t>{0U, 0U, 1U},
        "Texture IDs must follow first material/slot appearance and reuse path aliases.");

    expect(first.textures.size() == second.textures.size(),
           "Repeated packing must preserve the texture count.");
    expect(first.textureMipLevels.size() == second.textureMipLevels.size(),
           "Repeated packing must preserve the mip descriptor count.");
    for (std::size_t index = 0U;
         index < std::min(first.textures.size(), second.textures.size());
         ++index) {
        expect(
            sameTexture(first.textures[index], second.textures[index]),
            "Repeated packing must preserve stable texture descriptors.");
    }
    for (std::size_t index = 0U;
         index < std::min(
                     first.textureMipLevels.size(),
                     second.textureMipLevels.size());
         ++index) {
        expect(
            sameMip(first.textureMipLevels[index], second.textureMipLevels[index]),
            "Repeated packing must preserve stable mip descriptors.");
    }
    expect(
        first.textureTexelsRgba8 == second.textureTexelsRgba8,
        "Repeated packing must preserve every encoded RGBA8 texel.");
    expect(first.materials.size() == second.materials.size(),
           "Repeated packing must preserve the material count.");
    for (std::size_t index = 0U;
         index < std::min(first.materials.size(), second.materials.size());
         ++index) {
        expect(
            first.materials[index].textureIds ==
                second.materials[index].textureIds,
            "Repeated packing must preserve material texture IDs.");
    }
}

void testImportedEnvironmentPacking() {
    TemporaryDirectory directory;
    createImportedFixture(directory.path());
    writeTinyHdr(directory.path() / "tiny.hdr");

    Model model{directory.path(), "scene.obj", "tiny.hdr"};
    const PackedSceneData scene = model.packScene();
    scene.validate();

    expect(
        scene.environment.width == 2U && scene.environment.height == 1U &&
            scene.environment.flags == kPackedEnvironmentHasImportance,
        "The imported 2x1 HDR environment must publish sampleable metadata.");
    expect(
        scene.environmentRows.size() == 1U &&
            scene.environmentTexels.size() == 2U,
        "The imported HDR arrays must exactly match their extent.");
    if (scene.environmentRows.size() == 1U &&
        scene.environmentTexels.size() == 2U) {
        const PackedEnvironmentRow& row = scene.environmentRows.front();
        expect(
            row.probability == 1.0F && row.cumulativeProbability == 1.0F &&
                row.solidAngle > 0.0F,
            "A single HDR row must have unit probability and a positive solid angle.");
        const PackedEnvironmentTexel& red = scene.environmentTexels[0];
        const PackedEnvironmentTexel& green = scene.environmentTexels[1];
        expect(
            red.radianceAndConditionalCdf[0] > 0.0F &&
                red.radianceAndConditionalCdf[1] == 0.0F &&
                green.radianceAndConditionalCdf[0] == 0.0F &&
                green.radianceAndConditionalCdf[1] > 0.0F,
            "Packed HDR texels must preserve linear red and green radiance.");
        expect(
            std::abs(red.radianceAndConditionalCdf[3] - (1.0F / 3.0F)) <
                    1.0e-5F &&
                green.radianceAndConditionalCdf[3] == 1.0F,
            "The HDR conditional CDF must follow luminance times solid angle.");
    }
}

}  // namespace

int main() {
    try {
        testLightingValidation();
        testTextureValidation();
        testImportedTexturePacking();
        testImportedEnvironmentPacking();
    } catch (const std::exception& error) {
        std::cerr << "FAILED: unexpected exception: " << error.what() << '\n';
        return 1;
    }

    if (failureCount != 0) {
        std::cerr << failureCount << " scene-data test(s) failed.\n";
        return 1;
    }
    std::cout << "All scene-data tests passed.\n";
    return 0;
}
