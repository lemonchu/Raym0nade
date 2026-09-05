#include <chrono>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

#include "raym0nade/gpu/vulkan_path_renderer.hpp"
#include "raym0nade/model.hpp"
#include "raym0nade/render.hpp"
#include "raym0nade/scene_data.hpp"

#ifndef RAYM0NADE_PATH_TRACE_SHADER
#error "RAYM0NADE_PATH_TRACE_SHADER must be defined by CMake."
#endif
#ifndef RAYM0NADE_PATH_TRACE_SHADER_NAME
#error "RAYM0NADE_PATH_TRACE_SHADER_NAME must be defined by CMake."
#endif

namespace {

using Clock = std::chrono::steady_clock;
using namespace raym0nade;
using namespace raym0nade::gpu;

constexpr std::string_view kUnsupportedDeviceMessage =
    "No AMD Vulkan device satisfies the Ray Query backend requirements.";

struct Recipe {
    std::filesystem::path modelDirectory;
    std::filesystem::path modelFilename;
    std::filesystem::path skyFilename;
    RenderSettings settings;
};

struct Options {
    std::filesystem::path recipePath;
    std::filesystem::path shaderPath;
    std::optional<std::pair<int, int>> resolution;
    std::optional<int> samplesPerPixel;
    std::optional<std::uint32_t> seed;
    std::optional<float> exposure;
    std::optional<float> directLightProbability;
    std::optional<std::filesystem::path> outputPrefix;
    std::optional<std::size_t> expectedFaceCount;
    std::uint32_t tileWidth{128U};
    std::uint32_t tileHeight{128U};
    std::uint32_t samplesPerBatch{8U};
    bool requestValidation{false};
    bool exportAllPasses{false};
};

[[nodiscard]] double elapsedSeconds(const Clock::time_point begin) noexcept {
    return std::chrono::duration<double>(Clock::now() - begin).count();
}

[[noreturn]] void usage(const int exitCode) {
    std::ostream& stream = exitCode == 0 ? std::cout : std::cerr;
    stream
        << "Usage: raym0nade_gpu_render --recipe FILE [options]\n"
           "\n"
           "Reads a single-model, single-settings Raym0nade console recipe and renders it with Vulkan.\n"
           "\n"
           "Options:\n"
           "  --resolution WIDTHxHEIGHT  Override resolution and preserve horizontal view.\n"
           "  --spp N                    Override samples per pixel.\n"
           "  --seed N                   Override the deterministic uint32 seed.\n"
           "  --exposure X               Override the linear exposure multiplier.\n"
           "  --direct-probability X     Override direct-light allocation in [0, 1].\n"
           "  --output-prefix PATH       Override the recipe output prefix.\n"
           "  --expect-faces N           Fail before packing if import topology differs.\n"
           "  --tile WIDTHxHEIGHT        GPU tile extent (default: 128x128).\n"
           "  --batch-spp N              Samples per GPU dispatch, 1..64 (default: 8).\n"
           "  --shader FILE              Override the generated path-tracing SPIR-V.\n"
           "  --validation               Request Vulkan validation and synchronization checks.\n"
           "  --all-passes               Export the full CPU-compatible pass set.\n"
           "  --help                     Show this help.\n"
           "\n"
           "Relative paths use the current working directory, not the recipe directory.\n"
           "Model and sky filenames are resolved inside the model directory.\n"
           "Run checked-in examples from the repository root.\n";
    std::exit(exitCode);
}

[[nodiscard]] std::uint64_t parseUnsigned(
    const std::string& text, const char* optionName, const std::uint64_t maximum) {
    if (text.empty()) {
        throw std::invalid_argument(std::string{optionName} + " requires an unsigned integer.");
    }
    if (text.front() == '-') {
        throw std::invalid_argument(std::string{optionName} + " must not be negative.");
    }
    for (const char character : text) {
        if (character < '0' || character > '9') {
            throw std::invalid_argument(
                std::string{optionName} + " requires an unsigned decimal integer.");
        }
    }
    std::size_t consumed = 0U;
    unsigned long long value = 0ULL;
    try {
        value = std::stoull(text, &consumed, 10);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string{optionName} + " requires an unsigned integer.");
    }
    if (consumed != text.size() || value > maximum) {
        throw std::invalid_argument(std::string{optionName} + " is outside its supported range.");
    }
    return static_cast<std::uint64_t>(value);
}

[[nodiscard]] float parseFloat(const std::string& text, const char* optionName) {
    std::size_t consumed = 0U;
    float value = 0.0F;
    try {
        value = std::stof(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string{optionName} + " requires a finite number.");
    }
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument(std::string{optionName} + " requires a finite number.");
    }
    return value;
}

[[nodiscard]] std::pair<std::uint32_t, std::uint32_t> parseExtent(
    const std::string& text, const char* optionName) {
    const std::size_t separator = text.find_first_of("xX");
    if (separator == std::string::npos || separator == 0U ||
        separator + 1U >= text.size() || text.find_first_of("xX", separator + 1U) != std::string::npos) {
        throw std::invalid_argument(std::string{optionName} + " must use WIDTHxHEIGHT syntax.");
    }
    const std::uint64_t width = parseUnsigned(
        text.substr(0U, separator), optionName, std::numeric_limits<std::uint32_t>::max());
    const std::uint64_t height = parseUnsigned(
        text.substr(separator + 1U), optionName, std::numeric_limits<std::uint32_t>::max());
    if (width == 0U || height == 0U) {
        throw std::invalid_argument(std::string{optionName} + " dimensions must be positive.");
    }
    return {
        static_cast<std::uint32_t>(width),
        static_cast<std::uint32_t>(height),
    };
}

[[nodiscard]] std::filesystem::path resolveFromWorkingDirectory(
    const std::filesystem::path& path,
    const std::filesystem::path& workingDirectory) {
    if (path.is_absolute()) {
        return path.lexically_normal();
    }
    return (workingDirectory / path).lexically_normal();
}

[[nodiscard]] std::filesystem::path defaultShaderPath(
    const char* program,
    const std::filesystem::path& workingDirectory) {
    const std::filesystem::path programPath = program != nullptr
                                                  ? std::filesystem::path{program}
                                                  : std::filesystem::path{};
    const std::filesystem::path executable =
        resolveFromWorkingDirectory(programPath, workingDirectory);
    const std::filesystem::path adjacent =
        executable.parent_path() / RAYM0NADE_PATH_TRACE_SHADER_NAME;
    std::error_code error;
    if (std::filesystem::is_regular_file(adjacent, error)) {
        return adjacent;
    }
    return std::filesystem::path{RAYM0NADE_PATH_TRACE_SHADER};
}

[[nodiscard]] Options parseOptions(
    const int argc,
    char** argv,
    const std::filesystem::path& workingDirectory) {
    Options options;
    options.shaderPath = defaultShaderPath(argc > 0 ? argv[0] : nullptr, workingDirectory);
    for (int index = 1; index < argc; ++index) {
        const std::string argument{argv[index]};
        const auto next = [&](const char* optionName) -> std::string {
            if (index + 1 >= argc) {
                throw std::invalid_argument(std::string{optionName} + " requires a value.");
            }
            return argv[++index];
        };
        if (argument == "--recipe") {
            options.recipePath =
                resolveFromWorkingDirectory(next("--recipe"), workingDirectory);
        } else if (argument == "--resolution") {
            const auto extent = parseExtent(next("--resolution"), "--resolution");
            if (extent.first > static_cast<std::uint32_t>(std::numeric_limits<int>::max()) ||
                extent.second > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
                throw std::invalid_argument("--resolution exceeds the RenderSettings ABI.");
            }
            options.resolution = std::pair<int, int>{
                static_cast<int>(extent.first), static_cast<int>(extent.second)};
        } else if (argument == "--spp") {
            options.samplesPerPixel = static_cast<int>(parseUnsigned(
                next("--spp"), "--spp", static_cast<std::uint64_t>(std::numeric_limits<int>::max())));
        } else if (argument == "--seed") {
            options.seed = static_cast<std::uint32_t>(parseUnsigned(
                next("--seed"), "--seed", std::numeric_limits<std::uint32_t>::max()));
        } else if (argument == "--exposure") {
            options.exposure = parseFloat(next("--exposure"), "--exposure");
        } else if (argument == "--direct-probability") {
            options.directLightProbability =
                parseFloat(next("--direct-probability"), "--direct-probability");
        } else if (argument == "--output-prefix") {
            options.outputPrefix = next("--output-prefix");
        } else if (argument == "--expect-faces") {
            options.expectedFaceCount = static_cast<std::size_t>(parseUnsigned(
                next("--expect-faces"),
                "--expect-faces",
                std::numeric_limits<std::size_t>::max()));
        } else if (argument == "--tile") {
            const auto extent = parseExtent(next("--tile"), "--tile");
            options.tileWidth = extent.first;
            options.tileHeight = extent.second;
        } else if (argument == "--batch-spp") {
            options.samplesPerBatch = static_cast<std::uint32_t>(parseUnsigned(
                next("--batch-spp"),
                "--batch-spp",
                64U));
        } else if (argument == "--shader") {
            options.shaderPath =
                resolveFromWorkingDirectory(next("--shader"), workingDirectory);
        } else if (argument == "--validation") {
            options.requestValidation = true;
        } else if (argument == "--all-passes") {
            options.exportAllPasses = true;
        } else if (argument == "--help" || argument == "-h") {
            usage(0);
        } else {
            throw std::invalid_argument("Unknown option: " + argument);
        }
    }
    if (options.recipePath.empty()) {
        throw std::invalid_argument("--recipe is required.");
    }
    if (options.samplesPerPixel.has_value() && *options.samplesPerPixel <= 0) {
        throw std::invalid_argument("--spp must be positive.");
    }
    if (options.samplesPerBatch == 0U) {
        throw std::invalid_argument("--batch-spp must be positive.");
    }
    return options;
}

void requireToken(std::istream& input, const std::string& expected, const char* context) {
    std::string actual;
    if (!(input >> actual) || actual != expected) {
        throw std::runtime_error(
            std::string{"Invalid recipe while reading "} + context +
            ": expected '" + expected + "'.");
    }
}

template <typename Value>
void requireValue(std::istream& input, Value& value, const char* description) {
    if (!(input >> value)) {
        throw std::runtime_error(std::string{"Invalid recipe value for "} + description + '.');
    }
}

void readVector(std::istream& input, vec3& value, const char* description) {
    requireValue(input, value.x, description);
    requireValue(input, value.y, description);
    requireValue(input, value.z, description);
}

[[nodiscard]] Recipe loadRecipe(
    const std::filesystem::path& filename,
    const std::filesystem::path& workingDirectory) {
    std::ifstream input{filename};
    if (!input) {
        throw std::runtime_error("Could not open render recipe: " + filename.string());
    }

    Recipe result;
    std::string modelId;
    std::string settingsId;
    requireToken(input, "create", "model command");
    requireToken(input, "model", "model command");
    requireValue(input, modelId, "model identifier");
    requireValue(input, result.modelDirectory, "model directory");
    requireValue(input, result.modelFilename, "model filename");
    requireValue(input, result.skyFilename, "sky filename");

    requireToken(input, "create", "settings command");
    std::string settingsCommand;
    requireValue(input, settingsCommand, "settings command type");
    if (settingsCommand != "settings" && settingsCommand != "args") {
        throw std::runtime_error(
            "Invalid recipe while reading settings command: expected 'settings' or 'args'.");
    }
    requireValue(input, settingsId, "settings identifier");
    readVector(input, result.settings.direction, "camera direction");
    readVector(input, result.settings.right, "camera right vector");
    readVector(input, result.settings.up, "camera up vector");
    float alongDirection = 0.0F;
    float alongRight = 0.0F;
    float alongUp = 0.0F;
    requireValue(input, alongDirection, "camera position direction coefficient");
    requireValue(input, alongRight, "camera position right coefficient");
    requireValue(input, alongUp, "camera position up coefficient");
    result.settings.position = alongDirection * result.settings.direction +
                               alongRight * result.settings.right +
                               alongUp * result.settings.up;
    requireValue(input, result.settings.pixelScale, "pixel scale");
    requireValue(input, result.settings.focusDistance, "focus distance");
    requireValue(input, result.settings.circleOfConfusion, "circle of confusion");
    requireValue(input, result.settings.exposure, "exposure");
    requireValue(input, result.settings.width, "width");
    requireValue(input, result.settings.height, "height");
    requireValue(input, result.settings.samplesPerPixel, "samples per pixel");
    requireValue(input, result.settings.threadCount, "CPU thread count");
    requireValue(input, result.settings.directLightProbability, "direct-light probability");
    requireValue(input, result.settings.outputPrefix, "output prefix");

    requireToken(input, "render", "render command");
    std::string renderedModelId;
    std::string renderedSettingsId;
    requireValue(input, renderedModelId, "render model identifier");
    requireValue(input, renderedSettingsId, "render settings identifier");
    if (renderedModelId != modelId || renderedSettingsId != settingsId) {
        throw std::runtime_error("The recipe render command does not reference its created objects.");
    }
    std::string trailing;
    if (input >> trailing) {
        if (trailing != "exit") {
            throw std::runtime_error("Unexpected trailing recipe token: " + trailing);
        }
        if (input >> trailing) {
            throw std::runtime_error("Unexpected content after the recipe exit command.");
        }
    }
    result.modelDirectory =
        resolveFromWorkingDirectory(result.modelDirectory, workingDirectory);
    result.modelFilename = result.modelFilename.is_absolute()
                               ? result.modelFilename.lexically_normal()
                               : (result.modelDirectory / result.modelFilename).lexically_normal();
    if (!result.skyFilename.empty() && result.skyFilename != "null") {
        result.skyFilename = result.skyFilename.is_absolute()
                                 ? result.skyFilename.lexically_normal()
                                 : (result.modelFilename.parent_path() / result.skyFilename)
                                       .lexically_normal();
    }
    result.settings.outputPrefix =
        resolveFromWorkingDirectory(result.settings.outputPrefix, workingDirectory);
    return result;
}

void applyOverrides(
    const Options& options,
    const std::filesystem::path& workingDirectory,
    RenderSettings& settings) {
    if (options.resolution.has_value()) {
        const int previousWidth = settings.width;
        if (previousWidth <= 0) {
            throw std::invalid_argument(
                "Cannot apply --resolution because the recipe width is not positive.");
        }
        const int newWidth = options.resolution->first;
        const float scale = static_cast<float>(newWidth) / static_cast<float>(previousWidth);
        settings.pixelScale /= scale;
        settings.circleOfConfusion *= scale;
        settings.width = newWidth;
        settings.height = options.resolution->second;
    }
    if (options.samplesPerPixel.has_value()) {
        settings.samplesPerPixel = *options.samplesPerPixel;
    }
    if (options.seed.has_value()) {
        settings.seed = *options.seed;
    }
    if (options.exposure.has_value()) {
        settings.exposure = *options.exposure;
    }
    if (options.directLightProbability.has_value()) {
        settings.directLightProbability = *options.directLightProbability;
    }
    if (options.outputPrefix.has_value()) {
        settings.outputPrefix =
            resolveFromWorkingDirectory(*options.outputPrefix, workingDirectory);
    }
    settings.validate();
}

[[nodiscard]] std::filesystem::path beautyFilename(
    const std::filesystem::path& prefix) {
    std::filesystem::path filename = prefix;
    filename += "(Filter_FXAA).png";
    return filename;
}

void exportBeauty(Film film, const RenderSettings& settings) {
    const std::filesystem::path filename = beautyFilename(settings.outputPrefix);
    const std::filesystem::path parent = filename.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
    film.spatialClamp();
    film.filter();
    film.postProcess(Film::full | Film::fxaaEnabled);
    film.save(filename);
    std::cout << "Beauty image: " << filename.string() << '\n';
}

void printValidation(const VulkanValidationReport& report) {
    std::cout << "Validation: requested=" << (report.requested ? "yes" : "no")
              << ", enabled=" << (report.enabled ? "yes" : "no")
              << ", synchronization="
              << (report.synchronizationValidationEnabled ? "yes" : "no")
              << ", errors=" << report.errorCount
              << ", warnings=" << report.warningCount << '\n';
    for (const std::string& message : report.messages) {
        std::cout << "  " << message << '\n';
    }
}

}  // namespace

int main(const int argc, char** argv) {
    try {
        const std::filesystem::path workingDirectory = std::filesystem::current_path();
        const Options options = parseOptions(argc, argv, workingDirectory);
        Recipe recipe = loadRecipe(options.recipePath, workingDirectory);
        applyOverrides(options, workingDirectory, recipe.settings);

        const Clock::time_point totalBegin = Clock::now();
        VulkanPathRenderOptions renderOptions;
        renderOptions.vulkan.requestValidation = options.requestValidation;
        renderOptions.tileWidth = options.tileWidth;
        renderOptions.tileHeight = options.tileHeight;
        renderOptions.samplesPerBatch = options.samplesPerBatch;

        std::unique_ptr<VulkanPathRenderer> renderer;
        {
            PackedSceneData scene;
            {
                const Clock::time_point importBegin = Clock::now();
                Model model{recipe.modelDirectory, recipe.modelFilename, recipe.skyFilename};
                const double importSeconds = elapsedSeconds(importBegin);
                std::cout << "Imported " << model.faceCount() << " faces in " << std::fixed
                          << std::setprecision(3) << importSeconds << " s.\n";
                if (options.expectedFaceCount.has_value() &&
                    model.faceCount() != *options.expectedFaceCount) {
                    throw std::runtime_error(
                        "Imported face count does not match --expect-faces: expected " +
                        std::to_string(*options.expectedFaceCount) + ", got " +
                        std::to_string(model.faceCount()) + '.');
                }

                const Clock::time_point packBegin = Clock::now();
                scene = model.packScene();
                const double packSeconds = elapsedSeconds(packBegin);
                std::cout << "Packed " << scene.vertices.size() << " vertices, "
                          << scene.triangleCount() << " triangles, "
                          << scene.materials.size() << " materials, "
                          << scene.textures.size() << " textures, and "
                          << scene.textureTexelsRgba8.size() << " RGBA8 texels in "
                          << packSeconds << " s.\n";
            }
            renderer = std::make_unique<VulkanPathRenderer>(
                scene, options.shaderPath, renderOptions);
        }
        std::cout << "GPU: " << renderer->deviceName() << '\n';
        const VulkanRayQuerySetupTimings& setup = renderer->setupTimings();
        std::cout << "GPU upload: " << setup.uploadMilliseconds << " ms; AS build: "
                  << setup.accelerationStructureBuildMilliseconds << " ms.\n";

        VulkanPathRenderResult result = renderer->render(recipe.settings);
        std::cout << "GPU render host wall-clock: "
                  << result.timings.hostRenderMilliseconds << " ms.\n";
        if (result.timings.gpuTimestampAvailable) {
            std::cout << "GPU dispatch timestamps: "
                      << result.timings.gpuDispatchMilliseconds << " ms across "
                      << result.timings.dispatchCount << " dispatches.\n";
        }
        std::cout << "Direct-light events: " << result.stats.directLightSamples << '\n';

        const VulkanValidationReport validation = renderer->validationReport();
        printValidation(validation);
        if (options.requestValidation &&
            (!validation.enabled || !validation.synchronizationValidationEnabled ||
             validation.errorCount != 0U || validation.warningCount != 0U)) {
            throw std::runtime_error("Requested Vulkan validation did not complete cleanly.");
        }

        if (options.exportAllPasses) {
            exportFilmToFiles(std::move(result.film), recipe.settings);
            std::cout << "Exported all Film passes with prefix: "
                      << recipe.settings.outputPrefix.string() << '\n';
        } else {
            exportBeauty(std::move(result.film), recipe.settings);
        }
        std::cout << "Total wall-clock: " << elapsedSeconds(totalBegin) << " s.\n";
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "Error: " << exception.what() << '\n';
        if (std::string_view{exception.what()}.find(kUnsupportedDeviceMessage) !=
            std::string_view::npos) {
            return 77;
        }
        return 1;
    }
}
