#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include "raym0nade/gpu/vulkan_capabilities.hpp"
#include "raym0nade/gpu/vulkan_path_renderer.hpp"
#include "raym0nade/model.hpp"
#include "raym0nade/render.hpp"
#include "raym0nade/scene_data.hpp"

#include "render_recipe.hpp"

#ifndef RAYM0NADE_PATH_TRACE_SHADER
#error "RAYM0NADE_PATH_TRACE_SHADER must be defined by CMake."
#endif
#ifndef RAYM0NADE_PATH_TRACE_SHADER_NAME
#error "RAYM0NADE_PATH_TRACE_SHADER_NAME must be defined by CMake."
#endif
#ifndef RAYM0NADE_BENCHMARK_BUILD_TYPE
#define RAYM0NADE_BENCHMARK_BUILD_TYPE "unknown"
#endif

namespace {

using Clock = std::chrono::steady_clock;
using namespace raym0nade;
using namespace raym0nade::gpu;
using raym0nade::cli::SingleRenderRecipe;
using raym0nade::cli::loadSingleRenderRecipe;
using raym0nade::cli::resolveFromWorkingDirectory;

constexpr std::string_view kUnsupportedDeviceMessage =
    "No AMD Vulkan device satisfies the Ray Query backend requirements.";
constexpr std::size_t kMaximumWarmups = 100U;
constexpr std::size_t kMaximumMeasurements = 1000U;

struct Options {
    std::filesystem::path recipePath;
    std::filesystem::path shaderPath;
    std::optional<std::filesystem::path> comparisonShaderPath;
    std::filesystem::path outputDirectory{"output/benchmarks/path"};
    std::optional<std::pair<int, int>> resolution;
    std::optional<int> samplesPerPixel;
    std::optional<int> cpuThreadCount;
    std::optional<std::uint32_t> seed;
    std::optional<float> exposure;
    std::optional<float> directLightProbability;
    std::optional<std::size_t> expectedFaceCount;
    std::size_t warmupCount{1U};
    std::size_t measurementCount{3U};
    std::uint32_t tileWidth{128U};
    std::uint32_t tileHeight{128U};
    std::uint32_t samplesPerBatch{64U};
    std::uint32_t gpuQueueCount{1U};
    std::optional<std::uint32_t> comparisonGpuQueueCount;
    bool gpuOnly{false};
    bool unifiedCandidateGeometry{false};
    std::optional<bool> comparisonUnifiedCandidateGeometry;
    bool comparisonFirst{false};
    bool requestValidation{false};
};

struct IterationTiming {
    double externalWallMilliseconds{0.0};
    double reportedHostMilliseconds{0.0};
    std::optional<double> gpuTimestampMilliseconds;
    std::uint64_t directLightSamples{0U};
    std::uint64_t dispatchCount{0U};
};

struct BenchmarkRun {
    std::vector<IterationTiming> warmups;
    std::vector<IterationTiming> measurements;
    std::optional<Film> finalFilm;
};

struct Distribution {
    double minimum{0.0};
    double median{0.0};
    double percentile95{0.0};
    double maximum{0.0};
};

struct SceneTopology {
    std::size_t faces{0U};
    std::size_t vertices{0U};
    std::size_t triangles{0U};
    std::size_t materials{0U};
    std::size_t textures{0U};
    std::size_t textureTexels{0U};
    std::uint64_t textureTexelBytes{0U};
    std::size_t areaLights{0U};
    std::size_t areaLightTriangles{0U};
    std::size_t environmentTexels{0U};
    std::size_t cutoutTriangles{0U};
};

struct StageTimings {
    double importMilliseconds{0.0};
    double cpuExportMilliseconds{0.0};
    double packMilliseconds{0.0};
};

struct LinearQualityMetrics {
    double maximumAbsolute{0.0};
    double meanAbsolute{0.0};
    double rootMeanSquare{0.0};
    double cpuMeanLuminance{0.0};
    double gpuMeanLuminance{0.0};
};

struct ShaderArtifactIdentity {
    std::uint64_t byteCount{0U};
    std::uint64_t fnv1a64{0U};
};

struct CompletedGpuRun {
    std::filesystem::path shaderPath;
    ShaderArtifactIdentity shaderIdentity;
    std::string deviceName;
    VulkanRayQuerySetupTimings setup;
    VulkanValidationReport validation;
    BenchmarkRun benchmark;
    double constructionMilliseconds{0.0};
    double exportMilliseconds{0.0};
    std::optional<LinearQualityMetrics> cpuQuality;
    std::optional<std::vector<vec3>> retainedLinearBeauty;
    std::uint32_t computeQueueCount{1U};
    bool unifiedCandidateGeometry{false};
};

[[nodiscard]] constexpr const char* candidateGeometryDescription(
    const bool unified) noexcept {
    return unified ? "unified legacy ordering"
                   : "partitioned opaque/cutout ordering";
}

[[nodiscard]] constexpr const char* candidateGeometryToken(
    const bool unified) noexcept {
    return unified ? "unified" : "partitioned";
}

[[nodiscard]] constexpr const char* gpuExecutionOrderToken(
    const bool comparisonFirst) noexcept {
    return comparisonFirst ? "comparison-primary"
                           : "primary-comparison";
}

[[nodiscard]] double elapsedMilliseconds(
    const Clock::time_point begin) noexcept {
    return std::chrono::duration<double, std::milli>(
               Clock::now() - begin)
        .count();
}

[[noreturn]] void usage(const int exitCode) {
    std::ostream& stream = exitCode == 0 ? std::cout : std::cerr;
    stream
        << "Usage: raym0nade_gpu_path_benchmark --recipe FILE [options]\n"
           "\n"
           "Benchmarks the CPU and AMD Vulkan beauty path renderers from one imported Model.\n"
           "Use an optimized Release build for performance conclusions.\n"
           "\n"
           "Workload overrides:\n"
           "  --resolution WIDTHxHEIGHT  Override resolution and preserve horizontal view.\n"
           "  --spp N                    Override samples per pixel.\n"
           "  --cpu-threads N            CPU workers; 0 selects hardware concurrency.\n"
           "  --seed N                   Override the deterministic uint32 seed.\n"
           "  --exposure X               Override the linear exposure multiplier.\n"
           "  --direct-probability X     Override direct-light allocation in [0, 1].\n"
           "  --expect-faces N           Fail if the one import has a different face count.\n"
           "\n"
           "Measurement controls:\n"
           "  --warmups N                 Untimed-distribution warm-ups (default: 1; maximum: 100).\n"
           "  --measurements N            Measured renders (default: 3; maximum: 1000).\n"
           "  --tile WIDTHxHEIGHT         GPU tile extent (default: 128x128).\n"
           "  --batch-spp N               Samples per GPU dispatch, 1..64 (default: 64).\n"
           "  --gpu-queues N              Primary Vulkan compute queues (default: 1).\n"
           "  --comparison-gpu-queues N   Comparison-arm queues; also enables same-shader A/B.\n"
           "  --output-dir PATH           Reports and final PNGs (default: output/benchmarks/path).\n"
           "  --shader FILE               Override the generated path-tracing SPIR-V.\n"
           "  --comparison-shader FILE    Benchmark another SPIR-V from the same packed scene.\n"
           "  --unified-candidate-geometry Use legacy geometry ordering for the primary arm.\n"
           "  --comparison-unified-candidate-geometry Override comparison with legacy ordering.\n"
           "  --comparison-first          Run comparison before primary to reverse A/B order.\n"
           "  --gpu-only                  Skip CPU renders for repeated shader profiling.\n"
           "  --validation                Request Vulkan validation; not for timing claims.\n"
           "  --help                      Show this help.\n"
           "\n"
           "Defaults are smoke settings; use at least 10 measurements for performance work.\n"
           "PNG processing and export are outside every timed render call.\n"
           "Relative paths use the current working directory. Run repository recipes from the root.\n";
    std::exit(exitCode);
}

[[nodiscard]] std::uint64_t parseUnsigned(
    const std::string& text,
    const char* optionName,
    const std::uint64_t maximum) {
    if (text.empty()) {
        throw std::invalid_argument(
            std::string{optionName} + " requires an unsigned integer.");
    }
    for (const char character : text) {
        if (character < '0' || character > '9') {
            throw std::invalid_argument(
                std::string{optionName} +
                " requires an unsigned decimal integer.");
        }
    }
    std::size_t consumed = 0U;
    unsigned long long value = 0ULL;
    try {
        value = std::stoull(text, &consumed, 10);
    } catch (const std::exception&) {
        throw std::invalid_argument(
            std::string{optionName} + " requires an unsigned integer.");
    }
    if (consumed != text.size() || value > maximum) {
        throw std::invalid_argument(
            std::string{optionName} + " is outside its supported range.");
    }
    return static_cast<std::uint64_t>(value);
}

[[nodiscard]] ShaderArtifactIdentity identifyShaderArtifact(
    const std::filesystem::path& path) {
    std::ifstream input{path, std::ios::binary};
    if (!input) {
        throw std::runtime_error(
            "Could not open shader artifact for identity: " + path.string());
    }

    constexpr std::uint64_t kFnvOffsetBasis = 14695981039346656037ULL;
    constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
    std::array<char, 64U * 1024U> buffer{};
    ShaderArtifactIdentity identity{0U, kFnvOffsetBasis};
    for (;;) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize readCount = input.gcount();
        if (readCount < 0 ||
            static_cast<std::uint64_t>(readCount) >
                std::numeric_limits<std::uint64_t>::max() - identity.byteCount) {
            throw std::overflow_error("Shader artifact byte count overflowed.");
        }
        identity.byteCount += static_cast<std::uint64_t>(readCount);
        for (std::streamsize index = 0; index < readCount; ++index) {
            identity.fnv1a64 ^=
                static_cast<unsigned char>(buffer[static_cast<std::size_t>(index)]);
            identity.fnv1a64 *= kFnvPrime;
        }
        if (input.eof()) {
            break;
        }
        if (!input) {
            throw std::runtime_error(
                "Could not read shader artifact for identity: " + path.string());
        }
    }
    return identity;
}

[[nodiscard]] std::string describeShaderArtifact(
    const ShaderArtifactIdentity& identity) {
    std::ostringstream description;
    description << identity.byteCount
                << " bytes; FNV-1a-64=0x"
                << std::hex << std::nouppercase << std::setfill('0')
                << std::setw(16) << identity.fnv1a64;
    return description.str();
}

[[nodiscard]] std::string shaderChecksumToken(
    const ShaderArtifactIdentity& identity) {
    std::ostringstream token;
    token << std::hex << std::nouppercase << std::setfill('0')
          << std::setw(16) << identity.fnv1a64;
    return token.str();
}

[[nodiscard]] float parseFloat(
    const std::string& text,
    const char* optionName) {
    std::size_t consumed = 0U;
    float value = 0.0F;
    try {
        value = std::stof(text, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument(
            std::string{optionName} + " requires a finite number.");
    }
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument(
            std::string{optionName} + " requires a finite number.");
    }
    return value;
}

[[nodiscard]] std::pair<std::uint32_t, std::uint32_t> parseExtent(
    const std::string& text,
    const char* optionName) {
    const std::size_t separator = text.find_first_of("xX");
    if (separator == std::string::npos || separator == 0U ||
        separator + 1U >= text.size() ||
        text.find_first_of("xX", separator + 1U) != std::string::npos) {
        throw std::invalid_argument(
            std::string{optionName} + " must use WIDTHxHEIGHT syntax.");
    }
    const std::uint64_t width = parseUnsigned(
        text.substr(0U, separator),
        optionName,
        std::numeric_limits<std::uint32_t>::max());
    const std::uint64_t height = parseUnsigned(
        text.substr(separator + 1U),
        optionName,
        std::numeric_limits<std::uint32_t>::max());
    if (width == 0U || height == 0U) {
        throw std::invalid_argument(
            std::string{optionName} + " dimensions must be positive.");
    }
    return {
        static_cast<std::uint32_t>(width),
        static_cast<std::uint32_t>(height),
    };
}

[[nodiscard]] std::filesystem::path defaultShaderPath(
    const char* program,
    const std::filesystem::path& workingDirectory) {
    const std::filesystem::path programPath =
        program != nullptr ? std::filesystem::path{program}
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
    options.outputDirectory = resolveFromWorkingDirectory(
        options.outputDirectory,
        workingDirectory);
    options.shaderPath =
        defaultShaderPath(argc > 0 ? argv[0] : nullptr, workingDirectory);

    for (int index = 1; index < argc; ++index) {
        const std::string argument{argv[index]};
        const auto next = [&](const char* optionName) -> std::string {
            if (index + 1 >= argc) {
                throw std::invalid_argument(
                    std::string{optionName} + " requires a value.");
            }
            return argv[++index];
        };

        if (argument == "--recipe") {
            options.recipePath = resolveFromWorkingDirectory(
                next("--recipe"),
                workingDirectory);
        } else if (argument == "--resolution") {
            const auto extent =
                parseExtent(next("--resolution"), "--resolution");
            if (extent.first >
                    static_cast<std::uint32_t>(
                        std::numeric_limits<int>::max()) ||
                extent.second >
                    static_cast<std::uint32_t>(
                        std::numeric_limits<int>::max())) {
                throw std::invalid_argument(
                    "--resolution exceeds the RenderSettings ABI.");
            }
            options.resolution = {
                static_cast<int>(extent.first),
                static_cast<int>(extent.second),
            };
        } else if (argument == "--spp") {
            options.samplesPerPixel = static_cast<int>(parseUnsigned(
                next("--spp"),
                "--spp",
                static_cast<std::uint64_t>(
                    std::numeric_limits<int>::max())));
        } else if (argument == "--cpu-threads") {
            options.cpuThreadCount = static_cast<int>(parseUnsigned(
                next("--cpu-threads"),
                "--cpu-threads",
                static_cast<std::uint64_t>(
                    std::numeric_limits<int>::max())));
        } else if (argument == "--seed") {
            options.seed = static_cast<std::uint32_t>(parseUnsigned(
                next("--seed"),
                "--seed",
                std::numeric_limits<std::uint32_t>::max()));
        } else if (argument == "--exposure") {
            options.exposure =
                parseFloat(next("--exposure"), "--exposure");
        } else if (argument == "--direct-probability") {
            options.directLightProbability = parseFloat(
                next("--direct-probability"),
                "--direct-probability");
        } else if (argument == "--expect-faces") {
            options.expectedFaceCount =
                static_cast<std::size_t>(parseUnsigned(
                    next("--expect-faces"),
                    "--expect-faces",
                    std::numeric_limits<std::size_t>::max()));
        } else if (argument == "--warmups") {
            options.warmupCount =
                static_cast<std::size_t>(parseUnsigned(
                    next("--warmups"),
                    "--warmups",
                    kMaximumWarmups));
        } else if (argument == "--measurements") {
            options.measurementCount =
                static_cast<std::size_t>(parseUnsigned(
                    next("--measurements"),
                    "--measurements",
                    kMaximumMeasurements));
        } else if (argument == "--tile") {
            const auto extent = parseExtent(next("--tile"), "--tile");
            options.tileWidth = extent.first;
            options.tileHeight = extent.second;
        } else if (argument == "--batch-spp") {
            options.samplesPerBatch =
                static_cast<std::uint32_t>(parseUnsigned(
                    next("--batch-spp"),
                    "--batch-spp",
                    64U));
        } else if (argument == "--gpu-queues") {
            options.gpuQueueCount =
                static_cast<std::uint32_t>(parseUnsigned(
                    next("--gpu-queues"),
                    "--gpu-queues",
                    std::numeric_limits<std::uint32_t>::max()));
        } else if (argument == "--comparison-gpu-queues") {
            options.comparisonGpuQueueCount =
                static_cast<std::uint32_t>(parseUnsigned(
                    next("--comparison-gpu-queues"),
                    "--comparison-gpu-queues",
                    std::numeric_limits<std::uint32_t>::max()));
        } else if (argument == "--output-dir") {
            options.outputDirectory = resolveFromWorkingDirectory(
                next("--output-dir"),
                workingDirectory);
        } else if (argument == "--shader") {
            options.shaderPath = resolveFromWorkingDirectory(
                next("--shader"),
                workingDirectory);
        } else if (argument == "--comparison-shader") {
            options.comparisonShaderPath =
                resolveFromWorkingDirectory(
                    next("--comparison-shader"),
                    workingDirectory);
        } else if (argument == "--unified-candidate-geometry") {
            options.unifiedCandidateGeometry = true;
        } else if (argument ==
                   "--comparison-unified-candidate-geometry") {
            options.comparisonUnifiedCandidateGeometry = true;
        } else if (argument == "--comparison-first") {
            options.comparisonFirst = true;
        } else if (argument == "--gpu-only") {
            options.gpuOnly = true;
        } else if (argument == "--validation") {
            options.requestValidation = true;
        } else if (argument == "--help" || argument == "-h") {
            usage(0);
        } else {
            throw std::invalid_argument(
                "Unknown option: " + argument);
        }
    }

    if (options.recipePath.empty()) {
        throw std::invalid_argument("--recipe is required.");
    }
    if (options.outputDirectory.empty()) {
        throw std::invalid_argument("--output-dir must not be empty.");
    }
    if (options.measurementCount == 0U) {
        throw std::invalid_argument(
            "--measurements must be positive.");
    }
    if (options.samplesPerPixel.has_value() &&
        *options.samplesPerPixel <= 0) {
        throw std::invalid_argument("--spp must be positive.");
    }
    if (options.samplesPerBatch == 0U) {
        throw std::invalid_argument(
            "--batch-spp must be positive.");
    }
    if (options.gpuQueueCount == 0U) {
        throw std::invalid_argument("--gpu-queues must be positive.");
    }
    if (options.comparisonGpuQueueCount.has_value() &&
        *options.comparisonGpuQueueCount == 0U) {
        throw std::invalid_argument(
            "--comparison-gpu-queues must be positive.");
    }
    const bool hasComparisonArm =
        options.comparisonShaderPath.has_value() ||
        options.comparisonGpuQueueCount.has_value();
    if (!hasComparisonArm &&
        options.comparisonUnifiedCandidateGeometry) {
        throw std::invalid_argument(
            "--comparison-unified-candidate-geometry requires a comparison arm.");
    }
    if (!hasComparisonArm && options.comparisonFirst) {
        throw std::invalid_argument(
            "--comparison-first requires a comparison arm.");
    }
    return options;
}

void applyOverrides(
    const Options& options,
    RenderSettings& settings) {
    if (options.resolution.has_value()) {
        if (settings.width <= 0) {
            throw std::invalid_argument(
                "Cannot apply --resolution because the recipe width is not positive.");
        }
        const int oldWidth = settings.width;
        const int newWidth = options.resolution->first;
        const float scale =
            static_cast<float>(newWidth) / static_cast<float>(oldWidth);
        settings.pixelScale /= scale;
        settings.circleOfConfusion *= scale;
        settings.width = newWidth;
        settings.height = options.resolution->second;
    }
    if (options.samplesPerPixel.has_value()) {
        settings.samplesPerPixel = *options.samplesPerPixel;
    }
    if (options.cpuThreadCount.has_value()) {
        settings.threadCount = *options.cpuThreadCount;
    }
    if (options.seed.has_value()) {
        settings.seed = *options.seed;
    }
    if (options.exposure.has_value()) {
        settings.exposure = *options.exposure;
    }
    if (options.directLightProbability.has_value()) {
        settings.directLightProbability =
            *options.directLightProbability;
    }
    settings.validate();
}

[[nodiscard]] bool hasSupportedDevice() {
    const std::vector<VulkanDeviceCapabilities> devices =
        enumerateVulkanDevices();
    return std::any_of(
        devices.begin(),
        devices.end(),
        [](const VulkanDeviceCapabilities& device) {
            return device.supportsRayQueryBackend();
        });
}

[[nodiscard]] Distribution summarize(
    const std::vector<double>& values) {
    if (values.empty()) {
        throw std::invalid_argument(
            "Cannot summarize an empty timing sample.");
    }
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    const std::size_t middle = sorted.size() / 2U;
    const double median =
        sorted.size() % 2U == 0U
            ? (sorted[middle - 1U] + sorted[middle]) * 0.5
            : sorted[middle];
    const std::size_t percentile95Index =
        sorted.size() - 1U - sorted.size() / 20U;
    return {
        sorted.front(),
        median,
        sorted[percentile95Index],
        sorted.back(),
    };
}

[[nodiscard]] std::vector<double> selectExternalWall(
    const std::vector<IterationTiming>& timings) {
    std::vector<double> values;
    values.reserve(timings.size());
    for (const IterationTiming& timing : timings) {
        values.push_back(timing.externalWallMilliseconds);
    }
    return values;
}

[[nodiscard]] std::vector<double> selectReportedHost(
    const std::vector<IterationTiming>& timings) {
    std::vector<double> values;
    values.reserve(timings.size());
    for (const IterationTiming& timing : timings) {
        values.push_back(timing.reportedHostMilliseconds);
    }
    return values;
}

[[nodiscard]] std::vector<double> selectGpuTimestamps(
    const std::vector<IterationTiming>& timings) {
    std::vector<double> values;
    values.reserve(timings.size());
    for (const IterationTiming& timing : timings) {
        if (timing.gpuTimestampMilliseconds.has_value()) {
            values.push_back(*timing.gpuTimestampMilliseconds);
        }
    }
    return values;
}

[[nodiscard]] bool finiteVector(const vec3& value) noexcept {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
           std::isfinite(value.z);
}

[[nodiscard]] vec3 linearBeautyAt(
    const Film& film,
    const std::size_t index) {
    const RadianceData& directDiffuse =
        film.directDiffuseRadiance[index];
    const RadianceData& directSpecular =
        film.directSpecularRadiance[index];
    const RadianceData& indirectDiffuse =
        film.indirectDiffuseRadiance[index];
    const RadianceData& indirectSpecular =
        film.indirectSpecularRadiance[index];
    if (!finiteVector(directDiffuse.radiance) ||
        !finiteVector(directSpecular.radiance) ||
        !finiteVector(indirectDiffuse.radiance) ||
        !finiteVector(indirectSpecular.radiance) ||
        !std::isfinite(directDiffuse.varianceAccumulator) ||
        directDiffuse.varianceAccumulator < 0.0F ||
        !std::isfinite(directSpecular.varianceAccumulator) ||
        directSpecular.varianceAccumulator < 0.0F ||
        !std::isfinite(indirectDiffuse.varianceAccumulator) ||
        indirectDiffuse.varianceAccumulator < 0.0F ||
        !std::isfinite(indirectSpecular.varianceAccumulator) ||
        indirectSpecular.varianceAccumulator < 0.0F) {
        throw std::runtime_error(
            "A benchmark Film contains invalid radiance data.");
    }
    const HitInfo& hit = film.gBuffer[index];
    const vec3 value =
        hit.baseColor *
            (directDiffuse.radiance + indirectDiffuse.radiance) +
        directSpecular.radiance + indirectSpecular.radiance +
        hit.emission * film.exposure;
    if (!finiteVector(value)) {
        throw std::runtime_error(
            "A benchmark Film contains non-finite linear beauty.");
    }
    return value;
}

[[nodiscard]] std::vector<vec3> extractLinearBeauty(
    const Film& film) {
    const std::size_t pixelCount =
        static_cast<std::size_t>(film.width()) *
        static_cast<std::size_t>(film.height());
    if (film.gBuffer.size() != pixelCount ||
        film.directDiffuseRadiance.size() != pixelCount ||
        film.directSpecularRadiance.size() != pixelCount ||
        film.indirectDiffuseRadiance.size() != pixelCount ||
        film.indirectSpecularRadiance.size() != pixelCount) {
        throw std::runtime_error(
            "A benchmark Film has inconsistent plane sizes.");
    }

    std::vector<vec3> result(pixelCount, vec3{0.0F});
    for (std::size_t index = 0U; index < pixelCount; ++index) {
        result[index] = linearBeautyAt(film, index);
    }
    return result;
}

[[nodiscard]] double linearLuminance(const vec3& value) noexcept {
    return 0.3 * static_cast<double>(value.x) +
           0.6 * static_cast<double>(value.y) +
           0.1 * static_cast<double>(value.z);
}

[[nodiscard]] LinearQualityMetrics compareLinearBeauty(
    const std::vector<vec3>& reference,
    const std::vector<vec3>& candidate) {
    if (reference.size() != candidate.size()) {
        throw std::runtime_error(
            "Linear-beauty comparison extents differ.");
    }
    if (reference.empty()) {
        throw std::runtime_error(
            "Linear-beauty comparison has no pixels.");
    }

    LinearQualityMetrics metrics;
    double absoluteSum = 0.0;
    double squaredSum = 0.0;
    double cpuLuminanceSum = 0.0;
    double gpuLuminanceSum = 0.0;
    for (std::size_t index = 0U;
         index < reference.size();
         ++index) {
        const vec3& candidateValue = candidate[index];
        cpuLuminanceSum += linearLuminance(reference[index]);
        gpuLuminanceSum += linearLuminance(candidateValue);
        for (int channel = 0; channel < 3; ++channel) {
            const double difference = std::abs(
                static_cast<double>(candidateValue[channel]) -
                static_cast<double>(reference[index][channel]));
            metrics.maximumAbsolute =
                std::max(metrics.maximumAbsolute, difference);
            absoluteSum += difference;
            squaredSum += difference * difference;
        }
    }
    const double componentCount =
        static_cast<double>(reference.size()) * 3.0;
    metrics.meanAbsolute = absoluteSum / componentCount;
    metrics.rootMeanSquare =
        std::sqrt(squaredSum / componentCount);
    metrics.cpuMeanLuminance =
        cpuLuminanceSum / static_cast<double>(reference.size());
    metrics.gpuMeanLuminance =
        gpuLuminanceSum / static_cast<double>(reference.size());
    return metrics;
}

[[nodiscard]] SceneTopology inspectTopology(
    const PackedSceneData& scene,
    const std::size_t faceCount) {
    scene.validate();
    if (scene.textureTexelsRgba8.size() >
        std::numeric_limits<std::uint64_t>::max() /
            sizeof(std::uint32_t)) {
        throw std::overflow_error(
            "Packed texture byte count overflowed.");
    }

    SceneTopology topology;
    topology.faces = faceCount;
    topology.vertices = scene.vertices.size();
    topology.triangles = scene.triangleCount();
    topology.materials = scene.materials.size();
    topology.textures = scene.textures.size();
    topology.textureTexels = scene.textureTexelsRgba8.size();
    topology.textureTexelBytes =
        static_cast<std::uint64_t>(
            scene.textureTexelsRgba8.size()) *
        sizeof(std::uint32_t);
    topology.areaLights = scene.areaLights.size();
    topology.areaLightTriangles =
        scene.areaLightTriangles.size();
    topology.environmentTexels =
        scene.environmentTexels.size();
    for (const std::uint32_t materialId :
         scene.triangleMaterialIds) {
        const PackedMaterial& material =
            scene.materials[materialId];
        if ((material.flagsAndReserved[0] &
             kPackedMaterialCutout) != 0U) {
            ++topology.cutoutTriangles;
        }
    }
    return topology;
}

[[nodiscard]] BenchmarkRun benchmarkCpu(
    const Model& model,
    const RenderSettings& settings,
    const Options& options) {
    BenchmarkRun run;
    run.warmups.reserve(options.warmupCount);
    run.measurements.reserve(options.measurementCount);

    for (std::size_t iteration = 0U;
         iteration < options.warmupCount;
         ++iteration) {
        const Clock::time_point begin = Clock::now();
        FilmRenderResult output =
            renderToFilm(model, settings);
        const double externalMilliseconds =
            elapsedMilliseconds(begin);
        run.warmups.push_back({
            externalMilliseconds,
            output.stats.renderSeconds * 1.0e3,
            std::nullopt,
            output.stats.directLightSamples,
            0U,
        });
    }

    for (std::size_t iteration = 0U;
         iteration < options.measurementCount;
         ++iteration) {
        const Clock::time_point begin = Clock::now();
        FilmRenderResult output = renderToFilm(model, settings);
        const double externalMilliseconds =
            elapsedMilliseconds(begin);
        run.measurements.push_back({
            externalMilliseconds,
            output.stats.renderSeconds * 1.0e3,
            std::nullopt,
            output.stats.directLightSamples,
            0U,
        });
        run.finalFilm.emplace(std::move(output.film));
    }
    return run;
}

[[nodiscard]] BenchmarkRun benchmarkGpu(
    VulkanPathRenderer& renderer,
    const RenderSettings& settings,
    const Options& options) {
    BenchmarkRun run;
    run.warmups.reserve(options.warmupCount);
    run.measurements.reserve(options.measurementCount);

    const auto invoke = [&]() {
        const Clock::time_point begin = Clock::now();
        VulkanPathRenderResult output = renderer.render(settings);
        const double externalMilliseconds =
            elapsedMilliseconds(begin);
        IterationTiming timing{
            externalMilliseconds,
            output.timings.hostRenderMilliseconds,
            std::nullopt,
            output.stats.directLightSamples,
            output.timings.dispatchCount,
        };
        if (output.timings.gpuTimestampAvailable) {
            timing.gpuTimestampMilliseconds =
                output.timings.gpuDispatchMilliseconds;
        }
        return std::pair<IterationTiming, VulkanPathRenderResult>{
            timing,
            std::move(output),
        };
    };

    for (std::size_t iteration = 0U;
         iteration < options.warmupCount;
         ++iteration) {
        auto output = invoke();
        run.warmups.push_back(output.first);
    }
    for (std::size_t iteration = 0U;
         iteration < options.measurementCount;
         ++iteration) {
        auto output = invoke();
        run.measurements.push_back(output.first);
        run.finalFilm.emplace(std::move(output.second.film));
    }
    return run;
}

void requireStableWork(
    const BenchmarkRun& run,
    const char* backend,
    const bool requireDispatches) {
    if (run.measurements.empty()) {
        throw std::runtime_error(
            std::string{backend} + " produced no timing measurements.");
    }
    const IterationTiming& reference = run.measurements.front();
    for (const IterationTiming& timing : run.measurements) {
        if (timing.directLightSamples !=
            reference.directLightSamples) {
            throw std::runtime_error(
                std::string{backend} +
                " repeated renders executed different direct-light work.");
        }
        if (requireDispatches &&
            timing.dispatchCount != reference.dispatchCount) {
            throw std::runtime_error(
                std::string{backend} +
                " repeated renders executed different dispatch counts.");
        }
        if (timing.gpuTimestampMilliseconds.has_value() !=
            reference.gpuTimestampMilliseconds.has_value()) {
            throw std::runtime_error(
                std::string{backend} +
                " changed timestamp availability between measurements.");
        }
    }
}

[[nodiscard]] double exportBeauty(
    std::optional<Film>& film,
    const std::filesystem::path& filename) {
    if (!film.has_value()) {
        throw std::runtime_error(
            "The benchmark has no final Film to export.");
    }
    Film output = std::move(*film);
    film.reset();

    const Clock::time_point begin = Clock::now();
    output.spatialClamp();
    output.filter();
    output.postProcess(Film::full | Film::fxaaEnabled);
    output.save(filename);
    return elapsedMilliseconds(begin);
}

[[nodiscard]] std::size_t writeIterationRows(
    std::ostream& output,
    const char* backend,
    const char* phase,
    const std::vector<IterationTiming>& timings) {
    for (std::size_t index = 0U; index < timings.size(); ++index) {
        const IterationTiming& timing = timings[index];
        output << backend << ',' << phase << ',' << index << ','
               << timing.externalWallMilliseconds << ','
               << timing.reportedHostMilliseconds << ',';
        if (timing.gpuTimestampMilliseconds.has_value()) {
            output << *timing.gpuTimestampMilliseconds;
        }
        output << ',' << timing.directLightSamples << ','
               << timing.dispatchCount << '\n';
    }
    return timings.size();
}

[[nodiscard]] std::size_t writeTimingsCsv(
    const std::filesystem::path& filename,
    const BenchmarkRun* cpu,
    const BenchmarkRun& gpu,
    const BenchmarkRun* comparisonGpu) {
    std::ofstream output{filename};
    if (!output) {
        throw std::runtime_error(
            "Could not create benchmark timing CSV: " +
            filename.string());
    }
    output
        << "backend,phase,iteration,external_wall_ms,reported_host_ms,"
           "gpu_timestamp_ms,direct_light_events,dispatch_count\n"
        << std::setprecision(12);
    std::size_t rowCount = 0U;
    if (cpu != nullptr) {
        rowCount += writeIterationRows(
            output,
            "cpu",
            "warmup",
            cpu->warmups);
        rowCount += writeIterationRows(
            output,
            "cpu",
            "measurement",
            cpu->measurements);
    }
    rowCount += writeIterationRows(
        output,
        "gpu",
        "warmup",
        gpu.warmups);
    rowCount += writeIterationRows(
        output,
        "gpu",
        "measurement",
        gpu.measurements);
    if (comparisonGpu != nullptr) {
        rowCount += writeIterationRows(
            output,
            "comparison-gpu",
            "warmup",
            comparisonGpu->warmups);
        rowCount += writeIterationRows(
            output,
            "comparison-gpu",
            "measurement",
            comparisonGpu->measurements);
    }
    output.close();
    if (!output) {
        throw std::runtime_error(
            "Could not write benchmark timing CSV: " +
            filename.string());
    }

    std::ifstream verification{filename};
    if (!verification) {
        throw std::runtime_error(
            "Could not reopen benchmark timing CSV: " +
            filename.string());
    }
    std::size_t lineCount = 0U;
    std::string line;
    while (std::getline(verification, line)) {
        ++lineCount;
    }
    if (verification.bad() || lineCount != rowCount + 1U) {
        throw std::runtime_error(
            "Benchmark timing CSV verification failed: " +
            filename.string());
    }
    return rowCount;
}

void writeDistribution(
    std::ostream& output,
    const char* label,
    const Distribution& distribution) {
    output << "  " << label << ": median "
           << distribution.median << " ms (p95 "
           << distribution.percentile95 << ", min "
           << distribution.minimum << ", max "
           << distribution.maximum << ")\n";
}

[[nodiscard]] CompletedGpuRun executeGpuBenchmark(
    std::optional<PackedSceneData>& packedScene,
    const std::filesystem::path& shaderPath,
    const RenderSettings& settings,
    const Options& options,
    const std::filesystem::path& outputFilename,
    const std::vector<vec3>* cpuLinearBeauty,
    const bool unifiedCandidateGeometry,
    const std::uint32_t computeQueueCount,
    const bool retainLinearBeauty,
    const bool releasePackedAfterConstruction,
    const char* label) {
    if (!packedScene.has_value()) {
        throw std::logic_error(
            "Packed scene was released before GPU construction.");
    }

    const ShaderArtifactIdentity shaderIdentity =
        identifyShaderArtifact(shaderPath);

    VulkanPathRenderOptions renderOptions;
    renderOptions.vulkan.requestValidation =
        options.requestValidation;
    renderOptions.vulkan.forceUnifiedCandidateGeometry =
        unifiedCandidateGeometry;
    renderOptions.vulkan.computeQueueCount = computeQueueCount;
    renderOptions.tileWidth = options.tileWidth;
    renderOptions.tileHeight = options.tileHeight;
    renderOptions.samplesPerBatch =
        options.samplesPerBatch;

    const Clock::time_point constructionBegin = Clock::now();
    auto renderer = std::make_unique<VulkanPathRenderer>(
        *packedScene,
        shaderPath,
        renderOptions);
    const double constructionMilliseconds =
        elapsedMilliseconds(constructionBegin);
    if (releasePackedAfterConstruction) {
        packedScene.reset();
    }

    std::cout << label << " device: "
              << renderer->deviceName() << '\n'
              << label << " shader: "
              << shaderPath.string() << '\n'
              << label << " shader artifact: "
              << describeShaderArtifact(shaderIdentity) << '\n'
              << label << " candidate geometry: "
              << candidateGeometryDescription(
                     unifiedCandidateGeometry)
              << '\n'
              << label << " compute queues: "
              << renderer->computeQueueCount() << '\n'
              << label << " warm-ups: "
              << options.warmupCount
              << "; measurements: "
              << options.measurementCount << '\n';

    BenchmarkRun benchmark =
        benchmarkGpu(*renderer, settings, options);
    requireStableWork(benchmark, label, true);

    const bool needLinearBeauty =
        cpuLinearBeauty != nullptr || retainLinearBeauty;
    std::optional<std::vector<vec3>> linearBeauty;
    if (needLinearBeauty) {
        linearBeauty.emplace(
            extractLinearBeauty(*benchmark.finalFilm));
    }

    std::optional<LinearQualityMetrics> cpuQuality;
    if (cpuLinearBeauty != nullptr) {
        cpuQuality.emplace(compareLinearBeauty(
            *cpuLinearBeauty,
            *linearBeauty));
    }
    const VulkanRayQuerySetupTimings setup =
        renderer->setupTimings();
    const VulkanValidationReport validation =
        renderer->validationReport();
    if (options.requestValidation &&
        (!validation.enabled ||
         !validation.synchronizationValidationEnabled ||
         validation.errorCount != 0U ||
         validation.warningCount != 0U)) {
        throw std::runtime_error(
            "Requested Vulkan validation did not complete cleanly.");
    }

    const double exportMilliseconds =
        exportBeauty(benchmark.finalFilm, outputFilename);
    if (!retainLinearBeauty) {
        linearBeauty.reset();
    }

    return {
        shaderPath,
        shaderIdentity,
        renderer->deviceName(),
        setup,
        validation,
        std::move(benchmark),
        constructionMilliseconds,
        exportMilliseconds,
        cpuQuality,
        std::move(linearBeauty),
        renderer->computeQueueCount(),
        unifiedCandidateGeometry,
    };
}

void writeGpuDistributions(
    std::ostream& output,
    const char* heading,
    const CompletedGpuRun& gpu) {
    output << heading << "\n";
    writeDistribution(
        output,
        "external render-call wall clock",
        summarize(selectExternalWall(
            gpu.benchmark.measurements)));
    writeDistribution(
        output,
        "renderer-reported host time",
        summarize(selectReportedHost(
            gpu.benchmark.measurements)));
    const std::vector<double> timestamps =
        selectGpuTimestamps(gpu.benchmark.measurements);
    if (!timestamps.empty()) {
        writeDistribution(
            output,
            gpu.computeQueueCount == 1U
                ? "GPU dispatch timestamps"
                : "aggregate GPU queue-busy timestamps",
            summarize(timestamps));
    } else {
        output << "  GPU dispatch timestamps: unavailable\n";
    }
}

void writeQuality(
    std::ostream& output,
    const char* heading,
    const LinearQualityMetrics& metrics) {
    const std::ios::fmtflags savedFlags = output.flags();
    const std::streamsize savedPrecision = output.precision();
    output << std::scientific << std::setprecision(9)
           << heading << "\n"
           << "  Maximum absolute component error: "
           << metrics.maximumAbsolute << '\n'
           << "  Mean absolute component error: "
           << metrics.meanAbsolute << '\n'
           << "  Component RMSE: "
           << metrics.rootMeanSquare << '\n'
           << "  Reference mean luminance: "
           << metrics.cpuMeanLuminance << '\n'
           << "  Candidate mean luminance: "
           << metrics.gpuMeanLuminance << '\n';
    if (metrics.cpuMeanLuminance > 0.0) {
        output << "  Candidate/reference mean luminance: "
               << metrics.gpuMeanLuminance /
                      metrics.cpuMeanLuminance
               << '\n';
    }
    output.flags(savedFlags);
    output.precision(savedPrecision);
}

[[nodiscard]] std::string makeSummary(
    const Options& options,
    const SingleRenderRecipe& recipe,
    const RenderSettings& settings,
    const SceneTopology& topology,
    const StageTimings& stages,
    const BenchmarkRun* cpu,
    const CompletedGpuRun& gpu,
    const CompletedGpuRun* comparisonGpu,
    const LinearQualityMetrics* gpuComparisonQuality,
    const std::size_t timingRowCount,
    const double totalMilliseconds) {
    const Distribution gpuExternal = summarize(
        selectExternalWall(gpu.benchmark.measurements));
    const double pixelSamples =
        static_cast<double>(settings.width) *
        static_cast<double>(settings.height) *
        static_cast<double>(settings.samplesPerPixel);
    const double cutoutRatio =
        topology.triangles == 0U
            ? 0.0
            : static_cast<double>(topology.cutoutTriangles) /
                  static_cast<double>(topology.triangles);

    std::ostringstream output;
    output << std::fixed << std::setprecision(6)
           << "Raym0nade CPU/GPU beauty path benchmark\n\n"
           << "Build type: " << RAYM0NADE_BENCHMARK_BUILD_TYPE << '\n'
           << "Mode: "
           << (options.gpuOnly
                   ? "GPU-only shader profiling"
                   : "CPU/GPU same-import comparison")
           << '\n'
           << "Recipe: " << options.recipePath.string() << '\n'
           << "Model: " << recipe.modelFilename.string() << '\n'
           << "Environment: " << recipe.skyFilename.string() << '\n'
           << "Primary shader: " << gpu.shaderPath.string() << '\n'
           << "Primary shader artifact: "
           << describeShaderArtifact(gpu.shaderIdentity) << '\n'
           << "Primary candidate geometry: "
           << candidateGeometryDescription(
                  gpu.unifiedCandidateGeometry)
           << '\n'
           << "Primary compute queues: "
           << gpu.computeQueueCount << '\n';
    if (comparisonGpu != nullptr) {
        output << "Comparison shader: "
               << comparisonGpu->shaderPath.string() << '\n'
               << "Comparison shader artifact: "
               << describeShaderArtifact(comparisonGpu->shaderIdentity) << '\n'
               << "Comparison candidate geometry: "
               << candidateGeometryDescription(
                      comparisonGpu->unifiedCandidateGeometry)
               << '\n'
               << "Comparison compute queues: "
               << comparisonGpu->computeQueueCount << '\n'
               << "GPU execution order: "
               << gpuExecutionOrderToken(
                      options.comparisonFirst)
               << '\n';
    }

    output << "\nTopology from the one imported Model\n"
           << "  Faces: " << topology.faces << '\n'
           << "  Packed vertices: " << topology.vertices << '\n'
           << "  Packed triangles: " << topology.triangles << '\n'
           << "  Materials: " << topology.materials << '\n'
           << "  Textures: " << topology.textures << '\n'
           << "  Encoded RGBA8 texels: "
           << topology.textureTexels << '\n'
           << "  Encoded RGBA8 texel bytes: "
           << topology.textureTexelBytes << '\n'
           << "  Area lights: " << topology.areaLights << '\n'
           << "  Area-light triangles: "
           << topology.areaLightTriangles << '\n'
           << "  Environment texels: "
           << topology.environmentTexels << '\n'
           << "  Cutout triangles: "
           << topology.cutoutTriangles << " ("
           << cutoutRatio * 100.0 << "%)\n"
           << "\nFixed render settings\n"
           << "  Resolution: " << settings.width << 'x'
           << settings.height << '\n'
           << "  Samples per pixel: "
           << settings.samplesPerPixel << '\n'
           << "  Seed: " << settings.seed << '\n'
           << "  Exposure: " << settings.exposure << '\n'
           << "  Direct-light probability: "
           << settings.directLightProbability << '\n'
           << "  CPU thread setting / resolved: "
           << settings.threadCount << " / "
           << settings.resolvedThreadCount() << '\n'
           << "  GPU tile: " << options.tileWidth << 'x'
           << options.tileHeight << '\n'
           << "  GPU samples per dispatch: "
           << options.samplesPerBatch << '\n'
           << "  Warm-ups / measurements: "
           << options.warmupCount << " / "
           << options.measurementCount << '\n'
           << "  Measurement scope: "
           << (options.measurementCount < 10U
                   ? "smoke only"
                   : "performance-oriented sample count")
           << '\n'
           << "\nOne-time stages\n"
           << "  Model import and CPU BVH: "
           << stages.importMilliseconds << " ms\n";
    if (cpu != nullptr) {
        output << "  CPU final PNG export: "
               << stages.cpuExportMilliseconds << " ms\n";
    }
    output << "  Packed-scene conversion "
           << (cpu != nullptr ? "after CPU runs: " : "after import: ")
           << stages.packMilliseconds << " ms\n"
           << "  Source Model released before Vulkan construction: yes\n"
           << "  Primary Vulkan renderer construction: "
           << gpu.constructionMilliseconds << " ms\n"
           << "    Scene upload bucket: "
           << gpu.setup.uploadMilliseconds << " ms\n"
           << "    BLAS/TLAS build: "
           << gpu.setup.accelerationStructureBuildMilliseconds
           << " ms\n"
           << "  Primary GPU final PNG export: "
           << gpu.exportMilliseconds << " ms\n";
    if (comparisonGpu != nullptr) {
        output << "  Comparison Vulkan renderer construction: "
               << comparisonGpu->constructionMilliseconds
               << " ms\n"
               << "    Scene upload bucket: "
               << comparisonGpu->setup.uploadMilliseconds
               << " ms\n"
               << "    BLAS/TLAS build: "
               << comparisonGpu->setup
                      .accelerationStructureBuildMilliseconds
               << " ms\n"
               << "  Comparison GPU final PNG export: "
               << comparisonGpu->exportMilliseconds << " ms\n"
               << "  Packed host scene retained through both GPU measurement arms: yes\n"
               << "  Packed host scene released after both GPU arms returned: yes\n"
               << "  Two Vulkan renderers/device scenes alive concurrently: no\n";
    } else {
        output << "  Packed host scene released before primary GPU runs: yes\n";
    }

    if (cpu != nullptr) {
        output << "\nCPU measured distributions\n";
        const Distribution cpuExternal = summarize(
            selectExternalWall(cpu->measurements));
        writeDistribution(
            output,
            "external render-call wall clock",
            cpuExternal);
        writeDistribution(
            output,
            "renderer-reported kernel time",
            summarize(selectReportedHost(
                cpu->measurements)));
        output << "  Million pixel-samples/s at median: "
               << pixelSamples /
                      (cpuExternal.median * 1.0e3)
               << '\n'
               << "  Direct-light events per measured render: "
               << cpu->measurements.front().directLightSamples
               << '\n';
    }

    output << '\n';
    writeGpuDistributions(
        output,
        "Primary GPU measured distributions",
        gpu);
    output << "  Million pixel-samples/s at median: "
           << pixelSamples /
                  (gpuExternal.median * 1.0e3)
           << '\n'
           << "  Direct-light events per measured render: "
           << gpu.benchmark.measurements.front()
                  .directLightSamples
           << '\n'
           << "  Dispatches per measured render: "
           << gpu.benchmark.measurements.front().dispatchCount
           << '\n';

    if (cpu != nullptr) {
        const Distribution cpuExternal = summarize(
            selectExternalWall(cpu->measurements));
        output << "\nCPU/GPU host median speed ratio: "
               << cpuExternal.median / gpuExternal.median
               << "x\n"
               << "Values above 1 mean the primary GPU render call is faster.\n";
    }

    if (comparisonGpu != nullptr) {
        const Distribution comparisonExternal = summarize(
            selectExternalWall(
                comparisonGpu->benchmark.measurements));
        output << '\n';
        writeGpuDistributions(
            output,
            "Comparison GPU measured distributions",
            *comparisonGpu);
        output << "  Million pixel-samples/s at median: "
               << pixelSamples /
                      (comparisonExternal.median * 1.0e3)
               << '\n'
               << "  Direct-light events per measured render: "
               << comparisonGpu->benchmark.measurements.front()
                      .directLightSamples
               << '\n'
               << "  Dispatches per measured render: "
               << comparisonGpu->benchmark.measurements.front()
                      .dispatchCount
               << '\n'
               << "\nPrimary/comparison GPU host median ratio: "
               << gpuExternal.median /
                      comparisonExternal.median
               << "x\n"
               << "Values above 1 mean the comparison GPU arm is faster.\n";
    }

    if (gpu.cpuQuality.has_value()) {
        output << '\n';
        writeQuality(
            output,
            "CPU/primary-GPU linear beauty diagnostic",
            *gpu.cpuQuality);
    }
    if (comparisonGpu != nullptr &&
        comparisonGpu->cpuQuality.has_value()) {
        output << '\n';
        writeQuality(
            output,
            "CPU/comparison-GPU linear beauty diagnostic",
            *comparisonGpu->cpuQuality);
    }
    if (gpuComparisonQuality != nullptr) {
        output << '\n';
        writeQuality(
            output,
            "Primary/comparison-GPU linear beauty diagnostic",
            *gpuComparisonQuality);
    }

    output << "\nVulkan validation\n"
           << "  Requested / enabled / synchronization: "
           << (gpu.validation.requested ? "yes" : "no")
           << " / "
           << (gpu.validation.enabled ? "yes" : "no")
           << " / "
           << (gpu.validation.synchronizationValidationEnabled
                   ? "yes"
                   : "no")
           << '\n'
           << "  Errors / warnings: "
           << gpu.validation.errorCount << " / "
           << gpu.validation.warningCount << '\n';
    if (comparisonGpu != nullptr) {
        output << "  Comparison errors / warnings: "
               << comparisonGpu->validation.errorCount
               << " / "
               << comparisonGpu->validation.warningCount
               << '\n';
    }

    output << "\nTotal harness wall clock through PNG export: "
           << totalMilliseconds << " ms\n"
           << "\nTiming policy\n"
           << "  Warm-ups are recorded in timings.csv but excluded from distributions.\n"
           << "  Defaults of one warm-up and three measurements are smoke settings only.\n"
           << "  Use at least ten measurements before drawing performance conclusions.\n"
           << (cpu != nullptr
                   ? "  External render-call wall clock is the headline CPU/GPU comparison.\n"
                   : "  External render-call wall clock is the primary GPU host metric.\n")
           << "  GPU host time includes submission, waits, readback, and Film assembly.\n"
           << "  Readback is not reported as a separate timing bucket.\n"
           << "  GPU timestamps enclose device compute but exclude readback and Film assembly.\n"
           << "  Multi-queue timestamp totals are aggregate queue-busy time, not GPU wall clock.\n"
           << "  PNG processing and export are outside every render measurement.\n"
           << "  Linear beauty error is measured before spatial clamp, filter, gamma, and PNG encoding.\n"
           << "  CPU and GPU use different deterministic random streams, so this one-seed error\n"
           << "  is a stochastic diagnostic rather than a bias or golden-image threshold.\n"
           << "  Peak memory and rays per second are not measured; pixel-samples/s is reported.\n"
           << "  This harness alone does not implement the complete performance policy.\n";
    if (comparisonGpu != nullptr) {
        output
            << "  Fixed A/B order can be biased by driver warm-up, DVFS, and thermal state.\n"
            << "  Performance A/B requires separate AB and BA runs with matching topology.\n";
    }
    output << "  Verified timings.csv data rows: "
           << timingRowCount
           << "\nGenerated files\n";
    if (cpu != nullptr) {
        output << "  cpu-beauty.png\n";
    }
    output << "  gpu-beauty.png\n";
    if (comparisonGpu != nullptr) {
        output << "  comparison-gpu-beauty.png\n";
    }
    output << "  timings.csv\n"
           << "  summary.txt\n";
    return output.str();
}

}  // namespace

int main(const int argc, char** argv) {
    try {
        const std::filesystem::path workingDirectory =
            std::filesystem::current_path();
        const Options options =
            parseOptions(argc, argv, workingDirectory);
        if (std::string_view{RAYM0NADE_BENCHMARK_BUILD_TYPE} !=
            "Release") {
            std::cerr
                << "Warning: this benchmark was built as "
                << RAYM0NADE_BENCHMARK_BUILD_TYPE
                << "; use Release for performance conclusions.\n";
        }
        if (options.measurementCount < 10U) {
            std::cerr
                << "Note: fewer than ten measurements is a smoke run, "
                   "not a performance result.\n";
        }
        if (!hasSupportedDevice()) {
            std::cerr << "Skipped: " << kUnsupportedDeviceMessage
                      << '\n';
            return 77;
        }

        SingleRenderRecipe recipe = loadSingleRenderRecipe(
            options.recipePath,
            workingDirectory);
        applyOverrides(options, recipe.settings);
        std::filesystem::create_directories(
            options.outputDirectory);

        const Clock::time_point totalBegin = Clock::now();
        StageTimings stages;

        const Clock::time_point importBegin = Clock::now();
        auto model = std::make_unique<Model>(
            recipe.modelDirectory,
            recipe.modelFilename,
            recipe.skyFilename);
        stages.importMilliseconds =
            elapsedMilliseconds(importBegin);
        const std::size_t faceCount = model->faceCount();
        if (options.expectedFaceCount.has_value() &&
            faceCount != *options.expectedFaceCount) {
            throw std::runtime_error(
                "Imported face count does not match --expect-faces: expected " +
                std::to_string(*options.expectedFaceCount) +
                ", got " + std::to_string(faceCount) + '.');
        }
        std::cout << "Imported " << faceCount << " faces in "
                  << std::fixed << std::setprecision(3)
                  << stages.importMilliseconds << " ms.\n";

        std::optional<BenchmarkRun> cpu;
        std::optional<std::vector<vec3>> cpuLinearBeauty;
        if (!options.gpuOnly) {
            std::cout << "CPU warm-ups: " << options.warmupCount
                      << "; measurements: "
                      << options.measurementCount << '\n';
            cpu.emplace(
                benchmarkCpu(*model, recipe.settings, options));
            requireStableWork(*cpu, "CPU", false);
            cpuLinearBeauty.emplace(
                extractLinearBeauty(*cpu->finalFilm));
            stages.cpuExportMilliseconds = exportBeauty(
                cpu->finalFilm,
                options.outputDirectory / "cpu-beauty.png");
        }

        std::optional<PackedSceneData> packedScene;
        const Clock::time_point packBegin = Clock::now();
        packedScene.emplace(model->packScene());
        stages.packMilliseconds = elapsedMilliseconds(packBegin);
        const SceneTopology topology =
            inspectTopology(*packedScene, faceCount);
        model.reset();
        std::cout << "Packed " << topology.vertices
                  << " vertices, " << topology.triangles
                  << " triangles, " << topology.materials
                  << " materials, " << topology.textures
                  << " textures, and " << topology.textureTexels
                  << " RGBA8 texels in "
                  << stages.packMilliseconds << " ms.\n";

        const bool compareGpuArms =
            options.comparisonShaderPath.has_value() ||
            options.comparisonGpuQueueCount.has_value();
        const std::filesystem::path comparisonShaderPath =
            options.comparisonShaderPath.value_or(options.shaderPath);
        const std::uint32_t comparisonGpuQueueCount =
            options.comparisonGpuQueueCount.value_or(options.gpuQueueCount);
        const bool comparisonUnifiedCandidateGeometry =
            options.comparisonUnifiedCandidateGeometry.value_or(
                options.unifiedCandidateGeometry);
        std::optional<CompletedGpuRun> gpu;
        std::optional<CompletedGpuRun> comparisonGpu;
        const auto runPrimary = [&]() {
            return executeGpuBenchmark(
                packedScene,
                options.shaderPath,
                recipe.settings,
                options,
                options.outputDirectory / "gpu-beauty.png",
                cpuLinearBeauty.has_value()
                    ? &*cpuLinearBeauty
                    : nullptr,
                options.unifiedCandidateGeometry,
                options.gpuQueueCount,
                compareGpuArms,
                !compareGpuArms,
                "Primary GPU");
        };
        const auto runComparison = [&]() {
            return executeGpuBenchmark(
                packedScene,
                comparisonShaderPath,
                recipe.settings,
                options,
                options.outputDirectory /
                    "comparison-gpu-beauty.png",
                cpuLinearBeauty.has_value()
                    ? &*cpuLinearBeauty
                    : nullptr,
                comparisonUnifiedCandidateGeometry,
                comparisonGpuQueueCount,
                true,
                false,
                "Comparison GPU");
        };

        if (compareGpuArms && options.comparisonFirst) {
            comparisonGpu.emplace(runComparison());
            gpu.emplace(runPrimary());
        } else {
            gpu.emplace(runPrimary());
            if (compareGpuArms) {
                comparisonGpu.emplace(runComparison());
            }
        }

        std::optional<LinearQualityMetrics> gpuComparisonQuality;
        if (compareGpuArms) {
            if (!gpu->retainedLinearBeauty.has_value() ||
                !comparisonGpu->retainedLinearBeauty.has_value()) {
                throw std::logic_error(
                    "An A/B GPU linear beauty image was not retained.");
            }
            gpuComparisonQuality.emplace(compareLinearBeauty(
                *gpu->retainedLinearBeauty,
                *comparisonGpu->retainedLinearBeauty));
            gpu->retainedLinearBeauty.reset();
            comparisonGpu->retainedLinearBeauty.reset();
            packedScene.reset();
        }
        if (packedScene.has_value()) {
            throw std::logic_error(
                "Packed host scene outlived Vulkan construction.");
        }

        const double totalMilliseconds =
            elapsedMilliseconds(totalBegin);
        const std::filesystem::path csvPath =
            options.outputDirectory / "timings.csv";
        const std::size_t timingRowCount = writeTimingsCsv(
            csvPath,
            cpu.has_value() ? &*cpu : nullptr,
            gpu->benchmark,
            comparisonGpu.has_value()
                ? &comparisonGpu->benchmark
                : nullptr);

        const std::string summary = makeSummary(
            options,
            recipe,
            recipe.settings,
            topology,
            stages,
            cpu.has_value() ? &*cpu : nullptr,
            *gpu,
            comparisonGpu.has_value()
                ? &*comparisonGpu
                : nullptr,
            gpuComparisonQuality.has_value()
                ? &*gpuComparisonQuality
                : nullptr,
            timingRowCount,
            totalMilliseconds);
        const std::filesystem::path summaryPath =
            options.outputDirectory / "summary.txt";
        std::ofstream summaryFile{summaryPath};
        if (!summaryFile) {
            throw std::runtime_error(
                "Could not create benchmark summary: " +
                summaryPath.string());
        }
        summaryFile << summary;
        if (!summaryFile) {
            throw std::runtime_error(
                "Could not write benchmark summary: " +
                summaryPath.string());
        }

        std::cout << '\n' << summary
                  << "Reports written to: "
                  << options.outputDirectory.string() << '\n';
        if (comparisonGpu.has_value()) {
            const bool exactLinearMatch =
                gpuComparisonQuality->maximumAbsolute == 0.0 &&
                gpuComparisonQuality->meanAbsolute == 0.0 &&
                gpuComparisonQuality->rootMeanSquare == 0.0;
            std::cout
                << "Benchmark verification: timing-rows="
                << timingRowCount
                << " primary-layout="
                << candidateGeometryToken(
                       gpu->unifiedCandidateGeometry)
                << " comparison-layout="
                << candidateGeometryToken(
                       comparisonGpu->unifiedCandidateGeometry)
                << " linear-zero="
                << (exactLinearMatch ? "yes" : "no")
                << " order="
                << gpuExecutionOrderToken(
                       options.comparisonFirst)
                << " primary-queues="
                << gpu->computeQueueCount
                << " comparison-queues="
                << comparisonGpu->computeQueueCount
                << " primary-shader-bytes="
                << gpu->shaderIdentity.byteCount
                << " primary-shader-fnv1a64="
                << shaderChecksumToken(gpu->shaderIdentity)
                << " comparison-shader-bytes="
                << comparisonGpu->shaderIdentity.byteCount
                << " comparison-shader-fnv1a64="
                << shaderChecksumToken(comparisonGpu->shaderIdentity)
                << '\n';
        }
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "Error: " << exception.what() << '\n';
        if (std::string_view{exception.what()}.find(
                kUnsupportedDeviceMessage) !=
            std::string_view::npos) {
            return 77;
        }
        return 1;
    }
}
