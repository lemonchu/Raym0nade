#include <algorithm>
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
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "raym0nade/gpu/vulkan_primary_renderer.hpp"
#include "raym0nade/image.hpp"
#include "raym0nade/model.hpp"
#include "raym0nade/render.hpp"
#include "raym0nade/scene_data.hpp"

#ifndef RAYM0NADE_BENCHMARK_SOURCE_DIR
#error "RAYM0NADE_BENCHMARK_SOURCE_DIR must be defined by CMake."
#endif

#ifndef RAYM0NADE_PRIMARY_AOV_SHADER
#error "RAYM0NADE_PRIMARY_AOV_SHADER must be defined by CMake."
#endif

namespace {

using Clock = std::chrono::steady_clock;
using namespace raym0nade;
using namespace raym0nade::gpu;

constexpr float kAbsoluteTolerance = 1.0e-5F;
constexpr float kRelativeTolerance = 5.0e-5F;
constexpr float kDifferenceDisplayScale = 10000.0F;
constexpr std::size_t kExpectedBistroFaceCount = 2'832'120U;

struct Options {
    std::uint32_t width{1920U};
    std::uint32_t height{1080U};
    std::size_t warmupCount{3U};
    std::size_t measurementCount{10U};
    int cpuThreadCount{0};
    bool gpuOnly{false};
    bool requestValidation{false};
    bool forceUnifiedCandidateGeometry{false};
    std::filesystem::path outputDirectory{"output/benchmarks/g3b-bistro-shape-normal"};
};

struct Distribution {
    double minimum{0.0};
    double median{0.0};
    double percentile95{0.0};
    double maximum{0.0};
};

struct CpuRun {
    LinearImage image;
    double coldMilliseconds{0.0};
    std::vector<double> wallMilliseconds;
};

struct GpuRun {
    LinearImage image;
    double coldWallMilliseconds{0.0};
    std::vector<double> wallMilliseconds;
    std::vector<double> dispatchAndReadbackMilliseconds;
    std::vector<double> gpuDispatchMilliseconds;
};

struct ErrorMetrics {
    double maximumAbsolute{0.0};
    double meanAbsolute{0.0};
    double rootMeanSquare{0.0};
    std::size_t mismatchedPixels{0U};
    std::size_t nonFinitePixels{0U};
};

[[nodiscard]] double elapsedMilliseconds(const Clock::time_point start) noexcept {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

[[nodiscard]] std::uint32_t parseDimension(
    const std::string& text, const char* optionName) {
    const bool isDecimal = !text.empty() &&
                           std::all_of(text.begin(), text.end(), [](const char character) {
                               return character >= '0' && character <= '9';
                           });
    if (!isDecimal) {
        throw std::invalid_argument(
            std::string{optionName} + " must be a positive uint32 value.");
    }
    std::size_t consumed = 0U;
    const unsigned long long value = std::stoull(text, &consumed);
    if (consumed != text.size() || value == 0ULL ||
        value > static_cast<unsigned long long>(std::numeric_limits<std::uint32_t>::max())) {
        throw std::invalid_argument(std::string{optionName} + " must be a positive uint32 value.");
    }
    return static_cast<std::uint32_t>(value);
}

[[nodiscard]] std::size_t parseCount(const std::string& text, const char* optionName) {
    const bool isDecimal = !text.empty() &&
                           std::all_of(text.begin(), text.end(), [](const char character) {
                               return character >= '0' && character <= '9';
                           });
    if (!isDecimal) {
        throw std::invalid_argument(std::string{optionName} + " must be positive.");
    }
    std::size_t consumed = 0U;
    const unsigned long long value = std::stoull(text, &consumed);
    if (consumed != text.size() || value == 0ULL ||
        value > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max())) {
        throw std::invalid_argument(std::string{optionName} + " must be positive.");
    }
    return static_cast<std::size_t>(value);
}

[[nodiscard]] int parseThreadCount(const std::string& text) {
    std::size_t consumed = 0U;
    const long value = std::stol(text, &consumed);
    if (consumed != text.size() || value < 0L ||
        value > static_cast<long>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("--cpu-threads must be zero (automatic) or positive.");
    }
    return static_cast<int>(value);
}

[[noreturn]] void printUsageAndExit() {
    std::cout
        << "Usage: raym0nade_gpu_primary_benchmark [options]\n"
           "  --width N          Output width (default: 1920).\n"
           "  --height N         Output height (default: 1080).\n"
           "  --warmups N        Warm-up renders per active backend (default: 3).\n"
           "  --measurements N   Measured renders per active backend (default: 10).\n"
           "  --cpu-threads N    Parallel CPU workers; ignored with --gpu-only.\n"
           "  --output-dir PATH  Generated image and report directory.\n"
           "  --gpu-only         Skip CPU renders and CPU/GPU comparison outputs.\n"
           "  --unified-candidate-geometry\n"
           "                     Use the legacy single candidate geometry for GPU A/B.\n"
           "  --validation       Enable Vulkan validation (not for performance results).\n";
    std::exit(0);
}

[[nodiscard]] Options parseOptions(const int argc, char** argv) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string argument{argv[index]};
        const auto nextValue = [&](const char* optionName) -> std::string {
            if (index + 1 >= argc) {
                throw std::invalid_argument(std::string{optionName} + " requires a value.");
            }
            ++index;
            return argv[index];
        };

        if (argument == "--width") {
            options.width = parseDimension(nextValue("--width"), "--width");
        } else if (argument == "--height") {
            options.height = parseDimension(nextValue("--height"), "--height");
        } else if (argument == "--warmups") {
            options.warmupCount = parseCount(nextValue("--warmups"), "--warmups");
        } else if (argument == "--measurements") {
            options.measurementCount =
                parseCount(nextValue("--measurements"), "--measurements");
        } else if (argument == "--cpu-threads") {
            options.cpuThreadCount = parseThreadCount(nextValue("--cpu-threads"));
        } else if (argument == "--output-dir") {
            options.outputDirectory = nextValue("--output-dir");
        } else if (argument == "--gpu-only") {
            options.gpuOnly = true;
        } else if (argument == "--validation") {
            options.requestValidation = true;
        } else if (argument == "--unified-candidate-geometry") {
            options.forceUnifiedCandidateGeometry = true;
        } else if (argument == "--help" || argument == "-h") {
            printUsageAndExit();
        } else {
            throw std::invalid_argument("Unknown option: " + argument);
        }
    }
    return options;
}

[[nodiscard]] Distribution summarize(const std::vector<double>& values) {
    if (values.empty()) {
        throw std::invalid_argument("Cannot summarize an empty timing sample.");
    }
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    const std::size_t middle = sorted.size() / 2U;
    const double median = sorted.size() % 2U == 0U
                              ? (sorted[middle - 1U] + sorted[middle]) * 0.5
                              : sorted[middle];
    const std::size_t percentile95Index =
        sorted.size() - 1U - sorted.size() / 20U;
    return Distribution{
        sorted.front(), median, sorted[percentile95Index], sorted.back()};
}

[[nodiscard]] PrimaryRenderRequest makeRequest(
    const std::uint32_t width, const std::uint32_t height) {
    PrimaryRenderRequest request;
    request.extent = ImageExtent{width, height};

    // This is the checked-in Bistro exterior view. The legacy file stores position
    // coefficients along the camera basis; this is their resolved world position.
    request.camera.position = vec3{-15.594267F, 5.50455F, 2.0211134F};
    request.camera.direction = vec3{0.96718F, -0.2F, -0.156768F};
    request.camera.right = vec3{0.16F, 0.0F, 0.987117F};
    request.camera.up = vec3{-0.1974234F, -0.9798F, 0.032F};
    request.camera.pixelScale = 0.001025F * 2048.0F / static_cast<float>(width);
    request.aov = PrimaryAov::ShapeNormal;
    request.validate();
    return request;
}

[[nodiscard]] std::pair<LinearImage, double> renderCpuOnce(
    const Model& model,
    const PrimaryRenderRequest& request,
    const CpuPrimaryRenderOptions& options) {
    const auto start = Clock::now();
    LinearImage image = renderPrimaryAovCpu(model, request, options);
    image.validate();
    const double milliseconds = elapsedMilliseconds(start);
    return {std::move(image), milliseconds};
}

[[nodiscard]] CpuRun benchmarkCpu(
    const Model& model,
    const PrimaryRenderRequest& request,
    const CpuPrimaryRenderOptions& renderOptions,
    const Options& options) {
    auto cold = renderCpuOnce(model, request, renderOptions);
    CpuRun run{std::move(cold.first), cold.second, {}};

    for (std::size_t iteration = 0U; iteration < options.warmupCount; ++iteration) {
        run.image = renderCpuOnce(model, request, renderOptions).first;
    }
    run.wallMilliseconds.reserve(options.measurementCount);
    for (std::size_t iteration = 0U; iteration < options.measurementCount; ++iteration) {
        auto measured = renderCpuOnce(model, request, renderOptions);
        run.image = std::move(measured.first);
        run.wallMilliseconds.push_back(measured.second);
    }
    return run;
}

[[nodiscard]] std::pair<VulkanPrimaryRenderResult, double> renderGpuOnce(
    VulkanPrimaryRenderer& renderer, const PrimaryRenderRequest& request) {
    const auto start = Clock::now();
    VulkanPrimaryRenderResult result = renderer.render(request);
    const double milliseconds = elapsedMilliseconds(start);
    return {std::move(result), milliseconds};
}

[[nodiscard]] GpuRun benchmarkGpu(
    VulkanPrimaryRenderer& renderer,
    const PrimaryRenderRequest& request,
    const Options& options) {
    auto cold = renderGpuOnce(renderer, request);
    GpuRun run;
    run.image = std::move(cold.first.image);
    run.coldWallMilliseconds = cold.second;

    for (std::size_t iteration = 0U; iteration < options.warmupCount; ++iteration) {
        run.image = renderGpuOnce(renderer, request).first.image;
    }

    run.wallMilliseconds.reserve(options.measurementCount);
    run.dispatchAndReadbackMilliseconds.reserve(options.measurementCount);
    run.gpuDispatchMilliseconds.reserve(options.measurementCount);
    for (std::size_t iteration = 0U; iteration < options.measurementCount; ++iteration) {
        auto measured = renderGpuOnce(renderer, request);
        run.wallMilliseconds.push_back(measured.second);
        run.dispatchAndReadbackMilliseconds.push_back(
            measured.first.timings.dispatchAndReadbackMilliseconds);
        if (measured.first.timings.gpuTimestampAvailable) {
            run.gpuDispatchMilliseconds.push_back(
                measured.first.timings.gpuDispatchMilliseconds);
        }
        run.image = std::move(measured.first.image);
    }
    return run;
}

[[nodiscard]] std::size_t countNonFinitePixels(const LinearImage& image) noexcept {
    return static_cast<std::size_t>(std::count_if(
        image.pixels.begin(), image.pixels.end(), [](const vec3& pixel) {
            return !isFinite(pixel);
        }));
}

void requireExactlyEqual(const LinearImage& left, const LinearImage& right) {
    left.validate();
    right.validate();
    if (left.extent.width != right.extent.width ||
        left.extent.height != right.extent.height) {
        throw std::runtime_error("Single-thread and parallel CPU image extents differ.");
    }
    for (std::size_t index = 0U; index < left.pixels.size(); ++index) {
        const vec3& a = left.pixels[index];
        const vec3& b = right.pixels[index];
        if (a.x != b.x || a.y != b.y || a.z != b.z) {
            throw std::runtime_error(
                "Single-thread and parallel CPU images differ at pixel " +
                std::to_string(index) + '.');
        }
    }
}

[[nodiscard]] ErrorMetrics compareImages(
    const LinearImage& reference, const LinearImage& candidate) {
    reference.validate();
    candidate.validate();
    if (reference.extent.width != candidate.extent.width ||
        reference.extent.height != candidate.extent.height) {
        throw std::runtime_error("CPU and GPU image extents differ.");
    }

    ErrorMetrics metrics;
    double absoluteSum = 0.0;
    double squaredSum = 0.0;
    const std::size_t componentCount = reference.pixels.size() * 3U;
    for (std::size_t index = 0U; index < reference.pixels.size(); ++index) {
        bool mismatched = false;
        if (!isFinite(reference.pixels[index]) || !isFinite(candidate.pixels[index])) {
            ++metrics.nonFinitePixels;
            continue;
        }
        for (int channel = 0; channel < 3; ++channel) {
            const double expected = static_cast<double>(reference.pixels[index][channel]);
            const double actual = static_cast<double>(candidate.pixels[index][channel]);
            const double difference = std::abs(actual - expected);
            metrics.maximumAbsolute = std::max(metrics.maximumAbsolute, difference);
            absoluteSum += difference;
            squaredSum += difference * difference;
            const double tolerance = static_cast<double>(kAbsoluteTolerance) +
                                     static_cast<double>(kRelativeTolerance) *
                                         std::max(std::abs(expected), std::abs(actual));
            mismatched = mismatched || difference > tolerance;
        }
        if (mismatched) {
            ++metrics.mismatchedPixels;
        }
    }
    metrics.meanAbsolute = absoluteSum / static_cast<double>(componentCount);
    metrics.rootMeanSquare = std::sqrt(squaredSum / static_cast<double>(componentCount));
    return metrics;
}

[[nodiscard]] LinearImage makeDifferenceImage(
    const LinearImage& left, const LinearImage& right) {
    LinearImage difference{left.extent, std::vector<vec3>(left.pixels.size(), vec3{0.0F})};
    for (std::size_t index = 0U; index < left.pixels.size(); ++index) {
        difference.pixels[index] = glm::abs(left.pixels[index] - right.pixels[index]);
    }
    return difference;
}

[[nodiscard]] LinearImage makeComparisonImage(
    const LinearImage& cpu, const LinearImage& gpu) {
    constexpr std::uint32_t separatorWidth = 8U;
    if (cpu.extent.width >
        (std::numeric_limits<std::uint32_t>::max() - separatorWidth) / 2U) {
        throw std::length_error("Comparison image width is too large.");
    }
    const std::uint32_t width = cpu.extent.width * 2U + separatorWidth;
    LinearImage comparison{
        ImageExtent{width, cpu.extent.height},
        std::vector<vec3>(
            static_cast<std::size_t>(width) * cpu.extent.height, vec3{0.04F}),
    };
    for (std::uint32_t y = 0U; y < cpu.extent.height; ++y) {
        for (std::uint32_t x = 0U; x < cpu.extent.width; ++x) {
            const std::size_t source = static_cast<std::size_t>(y) * cpu.extent.width + x;
            const std::size_t row = static_cast<std::size_t>(y) * width;
            comparison.pixels[row + x] = cpu.pixels[source];
            comparison.pixels[row + cpu.extent.width + separatorWidth + x] = gpu.pixels[source];
        }
    }
    return comparison;
}

void writeTimingsCsv(
    const std::filesystem::path& path,
    const CpuRun& singleCpu,
    const CpuRun& parallelCpu,
    const GpuRun& gpu) {
    std::ofstream file{path};
    if (!file) {
        throw std::runtime_error("Could not create timing CSV: " + path.string());
    }
    file << "sample,cpu_single_wall_ms,cpu_parallel_wall_ms,gpu_wall_ms,"
            "gpu_dispatch_readback_ms,gpu_dispatch_ms\n";
    file << std::setprecision(9);
    for (std::size_t index = 0U; index < singleCpu.wallMilliseconds.size(); ++index) {
        file << index << ',' << singleCpu.wallMilliseconds[index] << ','
             << parallelCpu.wallMilliseconds[index] << ',' << gpu.wallMilliseconds[index]
             << ',' << gpu.dispatchAndReadbackMilliseconds[index] << ',';
        if (index < gpu.gpuDispatchMilliseconds.size()) {
            file << gpu.gpuDispatchMilliseconds[index];
        }
        file << '\n';
    }
}

void writeGpuOnlyTimingsCsv(
    const std::filesystem::path& path, const GpuRun& gpu) {
    std::ofstream file{path};
    if (!file) {
        throw std::runtime_error("Could not create timing CSV: " + path.string());
    }
    file << "mode,sample,gpu_wall_ms,gpu_dispatch_readback_ms,gpu_dispatch_ms\n";
    file << std::setprecision(9);
    for (std::size_t index = 0U; index < gpu.wallMilliseconds.size(); ++index) {
        file << "gpu-only," << index << ',' << gpu.wallMilliseconds[index] << ','
             << gpu.dispatchAndReadbackMilliseconds[index] << ',';
        if (index < gpu.gpuDispatchMilliseconds.size()) {
            file << gpu.gpuDispatchMilliseconds[index];
        }
        file << '\n';
    }
}

void printDistribution(
    std::ostream& output, const char* label, const std::vector<double>& values) {
    const Distribution distribution = summarize(values);
    output << "  " << label << ": median " << distribution.median << " ms"
           << " (p95 " << distribution.percentile95 << ", min "
           << distribution.minimum << ", max " << distribution.maximum << ")\n";
}

[[nodiscard]] std::string makeGpuOnlySummary(
    const Options& options,
    const Model& model,
    const PackedSceneData& scene,
    const VulkanPrimaryRenderer& renderer,
    const double modelMilliseconds,
    const double packMilliseconds,
    const double rendererMilliseconds,
    const GpuRun& gpu,
    const std::size_t nonFinitePixelCount) {
    const VulkanRayQuerySetupTimings& setup = renderer.setupTimings();
    const VulkanValidationReport validation = renderer.validationReport();

    std::ostringstream output;
    output << std::fixed << std::setprecision(6)
           << "Raym0nade GPU-only primary-AOV diagnostic\n\n"
           << "Mode: GPU-only; CPU rendering: skipped\n"
           << "Workload: ShapeNormal (one primary ray)\n"
           << "Resolution: " << options.width << 'x' << options.height << '\n'
           << "Scene: BistroExterior.fbx\n"
           << "Triangles: " << model.faceCount() << " (packed: " << scene.triangleCount()
           << ", expected complete Bistro: " << kExpectedBistroFaceCount << ")\n"
           << "Alpha cutout traversal: packed diffuse alpha candidate confirmation enabled\n"
           << "GPU: " << renderer.deviceName() << '\n'
           << "Warmups / measurements: " << options.warmupCount << " / "
           << options.measurementCount << '\n'
           << "Vulkan validation requested / enabled / synchronization: "
           << (validation.requested ? "yes" : "no") << " / "
           << (validation.enabled ? "yes" : "no") << " / "
           << (validation.synchronizationValidationEnabled ? "yes" : "no") << "\n\n"
           << "One-time scene and GPU setup\n"
           << "  Model import and CPU BVH (scene preparation only): "
           << modelMilliseconds << " ms\n"
           << "  Packed-scene conversion: " << packMilliseconds << " ms\n"
           << "  Vulkan renderer construction: " << rendererMilliseconds << " ms\n"
           << "    Scene upload bucket: " << setup.uploadMilliseconds << " ms\n"
           << "    BLAS/TLAS build: " << setup.accelerationStructureBuildMilliseconds
           << " ms\n\n"
           << "Cold wall clock\n"
           << "  GPU complete render call: " << gpu.coldWallMilliseconds << " ms\n\n"
           << "Warm GPU wall clock\n";
    printDistribution(output, "GPU complete render call", gpu.wallMilliseconds);
    printDistribution(output, "GPU dispatch and readback", gpu.dispatchAndReadbackMilliseconds);
    if (!gpu.gpuDispatchMilliseconds.empty()) {
        printDistribution(output, "GPU timestamp dispatch", gpu.gpuDispatchMilliseconds);
    }
    output << "\nGPU linear-image validation\n"
           << "  Non-finite pixels: " << nonFinitePixelCount << "\n\n"
           << "Generated by this run: gpu-shape-normal.png, timings.csv, summary.txt\n\n"
           << "Caveat: this is a GPU-only geometry diagnostic, not a textured or path-traced\n"
           << "beauty render. Packed diffuse alpha is sampled only for cutout visibility.\n";
    return output.str();
}

[[nodiscard]] std::string makeSummary(
    const Options& options,
    const Model& model,
    const PackedSceneData& scene,
    const VulkanPrimaryRenderer& renderer,
    const double modelMilliseconds,
    const double packMilliseconds,
    const double rendererMilliseconds,
    const int resolvedCpuThreads,
    const CpuRun& singleCpu,
    const CpuRun& parallelCpu,
    const GpuRun& gpu,
    const ErrorMetrics& errors) {
    const Distribution single = summarize(singleCpu.wallMilliseconds);
    const Distribution parallel = summarize(parallelCpu.wallMilliseconds);
    const Distribution gpuWall = summarize(gpu.wallMilliseconds);
    const VulkanRayQuerySetupTimings& setup = renderer.setupTimings();
    const VulkanValidationReport validation = renderer.validationReport();

    std::ostringstream output;
    output << std::fixed << std::setprecision(6)
           << "Raym0nade G3b CPU/GPU wall-clock comparison\n\n"
           << "Workload: ShapeNormal (one primary ray)\n"
           << "Resolution: " << options.width << 'x' << options.height << '\n'
           << "Scene: BistroExterior.fbx\n"
           << "Triangles: " << model.faceCount() << " (packed: " << scene.triangleCount()
           << ", expected complete Bistro: " << kExpectedBistroFaceCount << ")\n"
           << "Alpha cutout traversal: CPU/GPU packed diffuse alpha semantics enabled\n"
           << "GPU: " << renderer.deviceName() << '\n'
           << "Warmups / measurements: " << options.warmupCount << " / "
           << options.measurementCount << '\n'
           << "Vulkan validation requested / enabled / synchronization: "
           << (validation.requested ? "yes" : "no") << " / "
           << (validation.enabled ? "yes" : "no") << " / "
           << (validation.synchronizationValidationEnabled ? "yes" : "no") << "\n\n"
           << "One-time setup\n"
           << "  Model import and CPU BVH: " << modelMilliseconds << " ms\n"
           << "  Packed-scene conversion: " << packMilliseconds << " ms\n"
           << "  Vulkan renderer construction: " << rendererMilliseconds << " ms\n"
           << "    Scene upload bucket: " << setup.uploadMilliseconds << " ms\n"
           << "    BLAS/TLAS build: " << setup.accelerationStructureBuildMilliseconds << " ms\n\n"
           << "Cold wall clock\n"
           << "  CPU single-thread: " << singleCpu.coldMilliseconds << " ms\n"
           << "  CPU " << resolvedCpuThreads << " threads: "
           << parallelCpu.coldMilliseconds << " ms\n"
           << "  GPU complete render call: " << gpu.coldWallMilliseconds << " ms\n\n"
           << "Warm wall clock\n";
    printDistribution(output, "CPU single-thread", singleCpu.wallMilliseconds);
    output << "  CPU " << resolvedCpuThreads << " threads: median " << parallel.median
           << " ms (p95 " << parallel.percentile95 << ", min " << parallel.minimum
           << ", max " << parallel.maximum << ")\n";
    printDistribution(output, "GPU complete render call", gpu.wallMilliseconds);
    printDistribution(output, "GPU dispatch and readback", gpu.dispatchAndReadbackMilliseconds);
    if (!gpu.gpuDispatchMilliseconds.empty()) {
        printDistribution(output, "GPU timestamp dispatch", gpu.gpuDispatchMilliseconds);
    }
    output << "\nObserved median render-call wall-clock ratios\n"
           << "  CPU single-thread / GPU: " << single.median / gpuWall.median << "x\n"
           << "  CPU " << resolvedCpuThreads << " threads / GPU: "
           << parallel.median / gpuWall.median << "x\n\n"
           << "CPU/GPU linear-image difference\n"
           << "  Maximum absolute error: " << errors.maximumAbsolute << '\n'
           << "  Mean absolute error: " << errors.meanAbsolute << '\n'
           << "  RMSE: " << errors.rootMeanSquare << '\n'
           << "  Pixels outside abs+relative tolerance: " << errors.mismatchedPixels << '\n'
           << "  Non-finite pixels: " << errors.nonFinitePixels << "\n\n"
           << "Caveat: this ShapeNormal workload validates geometry and alpha-cutout traversal;\n"
           << "it is not a textured or path-traced beauty-render comparison.\n";
    return output.str();
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        std::filesystem::create_directories(options.outputDirectory);

        const std::filesystem::path sourceDirectory{RAYM0NADE_BENCHMARK_SOURCE_DIR};
        const std::filesystem::path sceneDirectory =
            sourceDirectory / "model" / "Bistro_v5_2";
        const auto modelStart = Clock::now();
        auto model = std::make_unique<Model>(
            sceneDirectory, "BistroExterior.fbx", "null");
        const double modelMilliseconds = elapsedMilliseconds(modelStart);
        if (model->faceCount() != kExpectedBistroFaceCount) {
            std::cerr << "Warning: Assimp imported " << model->faceCount()
                      << " Bistro faces; the complete asset has "
                      << kExpectedBistroFaceCount
                      << ". Treat this as an informal wall-clock run only.\n";
        }

        const auto packStart = Clock::now();
        PackedSceneData scene = model->packScene();
        const double packMilliseconds = elapsedMilliseconds(packStart);

        VulkanRayQueryOptions gpuOptions;
        gpuOptions.requestValidation = options.requestValidation;
        gpuOptions.forceUnifiedCandidateGeometry =
            options.forceUnifiedCandidateGeometry;
        const auto rendererStart = Clock::now();
        auto renderer = std::make_unique<VulkanPrimaryRenderer>(
            scene, std::filesystem::path{RAYM0NADE_PRIMARY_AOV_SHADER}, gpuOptions);
        const double rendererMilliseconds = elapsedMilliseconds(rendererStart);

        const PrimaryRenderRequest request = makeRequest(options.width, options.height);
        if (options.gpuOnly) {
            std::cout << "Rendering GPU-only ShapeNormal diagnostic at "
                      << options.width << 'x' << options.height << " on "
                      << renderer->deviceName() << "...\n";
            const GpuRun gpu = benchmarkGpu(*renderer, request, options);
            const std::size_t nonFinitePixelCount = countNonFinitePixels(gpu.image);
            if (nonFinitePixelCount == 0U) {
                saveLinearImagePng(
                    gpu.image, options.outputDirectory / "gpu-shape-normal.png");
            }
            writeGpuOnlyTimingsCsv(
                options.outputDirectory / "timings.csv", gpu);
            const std::string summary = makeGpuOnlySummary(
                options,
                *model,
                scene,
                *renderer,
                modelMilliseconds,
                packMilliseconds,
                rendererMilliseconds,
                gpu,
                nonFinitePixelCount);
            std::cout << '\n' << summary;
            std::ofstream summaryFile{options.outputDirectory / "summary.txt"};
            if (!summaryFile) {
                throw std::runtime_error("Could not create benchmark summary file.");
            }
            summaryFile << summary;

            if (options.requestValidation) {
                const VulkanValidationReport validation = renderer->validationReport();
                std::cout << "Validation enabled / errors / warnings: "
                          << (validation.enabled ? "yes" : "no") << " / "
                          << validation.errorCount << " / " << validation.warningCount << '\n';
                if (!validation.enabled || validation.errorCount != 0U ||
                    validation.warningCount != 0U) {
                    return 2;
                }
            }
            return nonFinitePixelCount == 0U ? 0 : 3;
        }

        const unsigned int hardwareThreads = std::thread::hardware_concurrency();
        const int resolvedCpuThreads = options.cpuThreadCount == 0
                                           ? static_cast<int>(std::max(hardwareThreads, 1U))
                                           : options.cpuThreadCount;

        std::cout << "Benchmarking " << options.width << 'x' << options.height << " on "
                  << renderer->deviceName() << "...\n"
                  << "CPU single-thread...\n";
        const CpuRun singleCpu = benchmarkCpu(
            *model, request, CpuPrimaryRenderOptions{1}, options);
        std::cout << "CPU " << resolvedCpuThreads << " threads...\n";
        const CpuRun parallelCpu = benchmarkCpu(
            *model, request, CpuPrimaryRenderOptions{options.cpuThreadCount}, options);
        requireExactlyEqual(singleCpu.image, parallelCpu.image);

        std::cout << "GPU...\n";
        const GpuRun gpu = benchmarkGpu(*renderer, request, options);
        const ErrorMetrics errors = compareImages(singleCpu.image, gpu.image);

        const LinearImage difference = makeDifferenceImage(singleCpu.image, gpu.image);
        const LinearImage comparison = makeComparisonImage(singleCpu.image, gpu.image);
        saveLinearImagePng(
            singleCpu.image, options.outputDirectory / "cpu-shape-normal.png");
        saveLinearImagePng(
            gpu.image, options.outputDirectory / "gpu-shape-normal.png");
        saveLinearImagePng(
            difference,
            options.outputDirectory / "absolute-error-x10000.png",
            kDifferenceDisplayScale);
        saveLinearImagePng(
            comparison, options.outputDirectory / "cpu-left-gpu-right.png");

        writeTimingsCsv(
            options.outputDirectory / "timings.csv", singleCpu, parallelCpu, gpu);
        const std::string summary = makeSummary(
            options,
            *model,
            scene,
            *renderer,
            modelMilliseconds,
            packMilliseconds,
            rendererMilliseconds,
            resolvedCpuThreads,
            singleCpu,
            parallelCpu,
            gpu,
            errors);
        std::cout << '\n' << summary;
        std::ofstream summaryFile{options.outputDirectory / "summary.txt"};
        if (!summaryFile) {
            throw std::runtime_error("Could not create benchmark summary file.");
        }
        summaryFile << summary;

        if (options.requestValidation) {
            const VulkanValidationReport validation = renderer->validationReport();
            std::cout << "Validation enabled / errors / warnings: "
                      << (validation.enabled ? "yes" : "no") << " / "
                      << validation.errorCount << " / " << validation.warningCount << '\n';
            if (!validation.enabled || validation.errorCount != 0U ||
                validation.warningCount != 0U) {
                return 2;
            }
        }
        return errors.nonFinitePixels == 0U ? 0 : 3;
    } catch (const std::exception& error) {
        std::cerr << "GPU primary benchmark failed: " << error.what() << '\n';
        return 1;
    }
}
