#include "raym0nade/gpu/vulkan_path_renderer.hpp"

#include <vulkan/vulkan.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "raym0nade/geometry.hpp"
#include "raym0nade/scene_data.hpp"
#include "vulkan_runtime.hpp"

namespace raym0nade::gpu {
namespace {

constexpr std::uint32_t kLocalSizeX = 8U;
constexpr std::uint32_t kLocalSizeY = 8U;
constexpr std::uint32_t kTimestampQueryCount = 2U;
constexpr std::uint32_t kStorageBufferBindingCount = 13U;
constexpr std::uint32_t kMaximumSamplesPerBatch = 64U;

struct alignas(16) PathRenderPushConstants {
    std::array<float, 4> cameraPositionAndPixelScale{};
    std::array<float, 4> cameraDirectionAndRayMinimum{};
    std::array<float, 4> cameraUpAndDirectProbability{};
    std::array<float, 4> cameraRightAndSampleWeight{};
    std::array<std::uint32_t, 4> imageExtentAndTileOrigin{};
    std::array<std::uint32_t, 4> tileExtentAndSampleRange{};
    std::array<std::uint32_t, 4> seedLightCountAndEnvironmentExtent{};
    std::array<std::uint32_t, 4> environmentAndSceneFlags{};
};

struct alignas(16) PathOutputPixel {
    std::array<float, 4> shapeNormalAndRoughness{};
    std::array<float, 4> surfaceNormalAndSpecular{};
    std::array<float, 4> positionAndMetallic{};
    std::array<float, 4> baseColorAndOpacity{};
    std::array<float, 4> emissionAndEta{};
    std::array<std::uint32_t, 4> materialEnteringHitAndDirectSampleCount{};
    std::array<float, 4> directDiffuse{};
    std::array<float, 4> directSpecular{};
    std::array<float, 4> indirectDiffuse{};
    std::array<float, 4> indirectSpecular{};
};

static_assert(sizeof(PathRenderPushConstants) == 128U);
static_assert(alignof(PathRenderPushConstants) == 16U);
static_assert(offsetof(PathRenderPushConstants, cameraPositionAndPixelScale) == 0U);
static_assert(offsetof(PathRenderPushConstants, cameraDirectionAndRayMinimum) == 16U);
static_assert(offsetof(PathRenderPushConstants, cameraUpAndDirectProbability) == 32U);
static_assert(offsetof(PathRenderPushConstants, cameraRightAndSampleWeight) == 48U);
static_assert(offsetof(PathRenderPushConstants, imageExtentAndTileOrigin) == 64U);
static_assert(offsetof(PathRenderPushConstants, tileExtentAndSampleRange) == 80U);
static_assert(offsetof(PathRenderPushConstants, seedLightCountAndEnvironmentExtent) == 96U);
static_assert(offsetof(PathRenderPushConstants, environmentAndSceneFlags) == 112U);
static_assert(std::is_standard_layout_v<PathRenderPushConstants>);
static_assert(std::is_trivially_copyable_v<PathRenderPushConstants>);

static_assert(sizeof(PathOutputPixel) == 160U);
static_assert(alignof(PathOutputPixel) == 16U);
static_assert(offsetof(PathOutputPixel, shapeNormalAndRoughness) == 0U);
static_assert(offsetof(PathOutputPixel, surfaceNormalAndSpecular) == 16U);
static_assert(offsetof(PathOutputPixel, positionAndMetallic) == 32U);
static_assert(offsetof(PathOutputPixel, baseColorAndOpacity) == 48U);
static_assert(offsetof(PathOutputPixel, emissionAndEta) == 64U);
static_assert(
    offsetof(PathOutputPixel, materialEnteringHitAndDirectSampleCount) == 80U);
static_assert(offsetof(PathOutputPixel, directDiffuse) == 96U);
static_assert(offsetof(PathOutputPixel, directSpecular) == 112U);
static_assert(offsetof(PathOutputPixel, indirectDiffuse) == 128U);
static_assert(offsetof(PathOutputPixel, indirectSpecular) == 144U);
static_assert(std::is_standard_layout_v<PathOutputPixel>);
static_assert(std::is_trivially_copyable_v<PathOutputPixel>);

[[nodiscard]] VulkanPathRenderOptions checkedOptions(
    VulkanPathRenderOptions options) {
    if (options.tileWidth == 0U || options.tileHeight == 0U) {
        throw std::invalid_argument("Vulkan path-render tile dimensions must be positive.");
    }
    if (options.samplesPerBatch == 0U) {
        throw std::invalid_argument(
            "Vulkan path-render samples per batch must be positive.");
    }
    if (options.samplesPerBatch > kMaximumSamplesPerBatch) {
        throw std::invalid_argument(
            "Vulkan path-render samples per batch must not exceed 64.");
    }
    return options;
}

[[nodiscard]] std::uint32_t checkedSceneCount(
    std::size_t count, const char* description) {
    if (count > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument(
            std::string{description} + " exceeds the GPU scene ABI limit.");
    }
    return static_cast<std::uint32_t>(count);
}

[[nodiscard]] std::uint32_t checkedMaterialCount(std::size_t count) {
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument(
            "The packed material count exceeds the Film material-ID limit.");
    }
    return static_cast<std::uint32_t>(count);
}

[[nodiscard]] std::uint32_t checkedGroupCount(
    std::uint32_t invocationCount,
    std::uint32_t localSize,
    std::uint32_t deviceLimit,
    const char* axis) {
    const std::uint64_t groups =
        (static_cast<std::uint64_t>(invocationCount) + localSize - 1U) / localSize;
    if (groups == 0U || groups > deviceLimit ||
        groups > std::numeric_limits<std::uint32_t>::max()) {
        throw std::invalid_argument(
            std::string{"The path render exceeds the Vulkan "} + axis +
            " dispatch limit.");
    }
    return static_cast<std::uint32_t>(groups);
}

[[nodiscard]] std::uint64_t timestampDelta(
    std::uint64_t begin, std::uint64_t end, std::uint32_t validBits) noexcept {
    if (validBits >= 64U) {
        return end - begin;
    }
    const std::uint64_t mask = (std::uint64_t{1U} << validBits) - 1U;
    return (end - begin) & mask;
}

[[nodiscard]] vec3 xyz(const std::array<float, 4>& value) noexcept {
    return vec3{value[0], value[1], value[2]};
}

void assignRadiance(
    const std::array<float, 4>& source,
    float exposure,
    RadianceData& destination) noexcept {
    destination.radiance = xyz(source);
    destination.varianceAccumulator = source[3];
    finalizeRadianceData(destination, exposure);
}

}  // namespace

class VulkanPathRenderer::Implementation {
public:
    Implementation(
        const PackedSceneData& scene,
        const std::filesystem::path& spirvPath,
        VulkanPathRenderOptions options)
        : options_(checkedOptions(options)),
          runtime_(scene, options_.vulkan),
          materialCount_(checkedMaterialCount(scene.materials.size())),
          areaLightCount_(checkedSceneCount(scene.areaLights.size(), "Area-light count")),
          environmentWidth_(scene.environment.width),
          environmentHeight_(scene.environment.height),
          environmentFlags_(scene.environment.flags),
          sceneVersion_(scene.formatVersion),
          queueStates_(runtime_.computeQueueCount()) {
        validateDeviceLimits();
        createPipeline(detail::readSpirvFile(spirvPath));
        createTimestampQueries();
    }

    [[nodiscard]] const std::string& deviceName() const noexcept {
        return runtime_.deviceName();
    }

    [[nodiscard]] std::uint32_t computeQueueCount() const noexcept {
        return runtime_.computeQueueCount();
    }

    [[nodiscard]] const VulkanRayQuerySetupTimings& setupTimings() const noexcept {
        return runtime_.setupTimings();
    }

    [[nodiscard]] VulkanValidationReport validationReport() const {
        return runtime_.validationReport();
    }

    [[nodiscard]] VulkanPathRenderResult render(const RenderSettings& settings) {
        settings.validate();

        const auto imageWidth = static_cast<std::uint32_t>(settings.width);
        const auto imageHeight = static_cast<std::uint32_t>(settings.height);
        const auto totalSamples =
            static_cast<std::uint32_t>(settings.samplesPerPixel);
        validateRenderArithmetic(imageWidth, imageHeight, totalSamples);

        const detail::VulkanClock::time_point renderBegin = detail::VulkanClock::now();
        Film film{settings.width, settings.height};
        film.exposure = settings.exposure;
        film.focusDistance = settings.focusDistance;
        film.circleOfConfusion = settings.circleOfConfusion;
        film.cameraPosition = settings.position;

        VulkanPathRenderTimings timings;
        timings.gpuTimestampAvailable = timestampQueriesAvailable_;
        timings.computeQueueCount = runtime_.computeQueueCount();

        std::lock_guard<std::mutex> lock{runtime_.operationMutex()};
        const std::uint32_t maximumTileWidth =
            std::min(options_.tileWidth, imageWidth);
        const std::uint32_t maximumTileHeight =
            std::min(options_.tileHeight, imageHeight);
        const std::size_t maximumTilePixels =
            static_cast<std::size_t>(maximumTileWidth) * maximumTileHeight;
        const VkDeviceSize maximumOutputBytes = detail::checkedVulkanByteSize(
            maximumTilePixels,
            sizeof(PathOutputPixel),
            "Vulkan path-render tile output");
        for (QueueState& queueState : queueStates_) {
            ensureOutputCapacity(
                queueState, maximumTilePixels, maximumOutputBytes);
            queueState.hostPixels.resize(maximumTilePixels);
        }

        const std::size_t tileColumns =
            (static_cast<std::size_t>(imageWidth) - 1U) / options_.tileWidth + 1U;
        const std::size_t tileRows =
            (static_cast<std::size_t>(imageHeight) - 1U) / options_.tileHeight + 1U;
        if (tileRows > std::numeric_limits<std::size_t>::max() / tileColumns) {
            throw std::overflow_error("The Vulkan path-render tile count overflowed.");
        }
        const std::size_t tileCount = tileColumns * tileRows;
        std::atomic<std::size_t> nextTile{0U};
        std::atomic<bool> stopWorkers{false};
        std::vector<WorkerResult> workerResults(queueStates_.size());

        const auto renderWorker = [&](const std::uint32_t queueIndex) noexcept {
            WorkerResult& worker = workerResults[queueIndex];
            try {
                while (!stopWorkers.load(std::memory_order_relaxed)) {
                    const std::size_t tileIndex =
                        nextTile.fetch_add(1U, std::memory_order_relaxed);
                    if (tileIndex >= tileCount) {
                        break;
                    }
                    const std::size_t tileColumn = tileIndex % tileColumns;
                    const std::size_t tileRow = tileIndex / tileColumns;
                    const auto tileX = static_cast<std::uint32_t>(
                        tileColumn * options_.tileWidth);
                    const auto tileY = static_cast<std::uint32_t>(
                        tileRow * options_.tileHeight);
                    const std::uint32_t tileWidth =
                        std::min(options_.tileWidth, imageWidth - tileX);
                    const std::uint32_t tileHeight =
                        std::min(options_.tileHeight, imageHeight - tileY);
                    renderTile(
                        queueIndex,
                        queueStates_[queueIndex],
                        settings,
                        imageWidth,
                        imageHeight,
                        tileX,
                        tileY,
                        tileWidth,
                        tileHeight,
                        totalSamples,
                        film,
                        worker.directLightSamples,
                        worker.gpuDispatchMilliseconds,
                        worker.dispatchCount);
                }
            } catch (...) {
                worker.exception = std::current_exception();
                stopWorkers.store(true, std::memory_order_relaxed);
            }
        };

        std::vector<std::thread> workers;
        workers.reserve(queueStates_.size() - 1U);
        try {
            for (std::uint32_t queueIndex = 1U;
                 queueIndex < queueStates_.size();
                 ++queueIndex) {
                workers.emplace_back(renderWorker, queueIndex);
            }
        } catch (...) {
            stopWorkers.store(true, std::memory_order_relaxed);
            for (std::thread& worker : workers) {
                worker.join();
            }
            throw;
        }
        renderWorker(0U);
        for (std::thread& worker : workers) {
            worker.join();
        }

        std::uint64_t directLightSamples = 0U;
        for (const WorkerResult& worker : workerResults) {
            if (worker.exception != nullptr) {
                std::rethrow_exception(worker.exception);
            }
            if (worker.directLightSamples >
                std::numeric_limits<std::uint64_t>::max() - directLightSamples) {
                throw std::overflow_error(
                    "The GPU direct-light sample count overflowed.");
            }
            directLightSamples += worker.directLightSamples;
            timings.gpuDispatchMilliseconds += worker.gpuDispatchMilliseconds;
            if (worker.dispatchCount >
                std::numeric_limits<std::uint64_t>::max() - timings.dispatchCount) {
                throw std::overflow_error("The GPU dispatch count overflowed.");
            }
            timings.dispatchCount += worker.dispatchCount;
        }

        timings.hostRenderMilliseconds = detail::elapsedMilliseconds(renderBegin);
        RenderStats stats;
        stats.renderSeconds = timings.hostRenderMilliseconds / 1.0e3;
        stats.totalSeconds = stats.renderSeconds;
        stats.directLightSamples = directLightSamples;
        return VulkanPathRenderResult{std::move(film), stats, timings};
    }

private:
    struct QueueState {
        std::size_t outputCapacity{0U};
        std::vector<PathOutputPixel> hostPixels;
        std::unique_ptr<detail::VulkanBuffer> outputBuffer;
        std::unique_ptr<detail::VulkanBuffer> readbackBuffer;
        VkDescriptorSet descriptorSet{VK_NULL_HANDLE};
        detail::UniqueVulkanHandle<VkQueryPool> timestampQueryPool;
    };

    struct WorkerResult {
        std::uint64_t directLightSamples{0U};
        double gpuDispatchMilliseconds{0.0};
        std::uint64_t dispatchCount{0U};
        std::exception_ptr exception;
    };

    void validateDeviceLimits() const {
        const VkPhysicalDeviceLimits& limits = runtime_.physicalProperties().limits;
        if (limits.maxComputeWorkGroupSize[0] < kLocalSizeX ||
            limits.maxComputeWorkGroupSize[1] < kLocalSizeY ||
            limits.maxComputeWorkGroupSize[2] < 1U ||
            limits.maxComputeWorkGroupInvocations < kLocalSizeX * kLocalSizeY) {
            throw std::runtime_error(
                "The Vulkan device cannot execute the path renderer's 8 x 8 workgroup.");
        }
        if (limits.maxComputeWorkGroupCount[2] < 1U) {
            throw std::runtime_error(
                "The Vulkan device cannot dispatch the path renderer in Z.");
        }
        if (limits.maxPushConstantsSize < sizeof(PathRenderPushConstants)) {
            throw std::runtime_error(
                "The Vulkan device does not provide 128 bytes of push-constant storage.");
        }
        if (limits.maxPerStageDescriptorStorageBuffers < kStorageBufferBindingCount ||
            limits.maxDescriptorSetStorageBuffers < kStorageBufferBindingCount) {
            throw std::runtime_error(
                "The Vulkan device does not provide thirteen storage-buffer descriptors.");
        }
        if (limits.maxPerStageResources < kStorageBufferBindingCount + 1U) {
            throw std::runtime_error(
                "The Vulkan device does not provide fourteen per-stage resources.");
        }

        const std::array<const detail::VulkanBuffer*, 12> sceneBuffers{{
            &runtime_.vertexBuffer(),
            &runtime_.indexBuffer(),
            &runtime_.triangleMaterialIdBuffer(),
            &runtime_.materialBuffer(),
            &runtime_.textureDescriptorBuffer(),
            &runtime_.textureMipBuffer(),
            &runtime_.textureTexelBuffer(),
            &runtime_.areaLightBuffer(),
            &runtime_.areaLightTriangleBuffer(),
            &runtime_.environmentRowBuffer(),
            &runtime_.environmentTexelBuffer(),
            &runtime_.primitiveRemapBuffer(),
        }};
        for (const detail::VulkanBuffer* buffer : sceneBuffers) {
            if (buffer->size() > limits.maxStorageBufferRange) {
                throw std::runtime_error(
                    "A packed-scene buffer exceeds the Vulkan storage-buffer range limit.");
            }
        }
    }

    void validateRenderArithmetic(
        std::uint32_t imageWidth,
        std::uint32_t imageHeight,
        std::uint32_t totalSamples) const {
        const VkPhysicalDeviceLimits& limits = runtime_.physicalProperties().limits;
        const std::uint32_t maximumTileWidth =
            std::min(options_.tileWidth, imageWidth);
        const std::uint32_t maximumTileHeight =
            std::min(options_.tileHeight, imageHeight);
        static_cast<void>(checkedGroupCount(
            maximumTileWidth,
            kLocalSizeX,
            limits.maxComputeWorkGroupCount[0],
            "X"));
        static_cast<void>(checkedGroupCount(
            maximumTileHeight,
            kLocalSizeY,
            limits.maxComputeWorkGroupCount[1],
            "Y"));

        const VkDeviceSize maximumOutputBytes = detail::checkedVulkanByteSize(
            static_cast<std::size_t>(maximumTileWidth) *
                static_cast<std::size_t>(maximumTileHeight),
            sizeof(PathOutputPixel),
            "Vulkan path-render tile output");
        if (maximumOutputBytes > limits.maxStorageBufferRange) {
            throw std::invalid_argument(
                "The configured path-render tile exceeds maxStorageBufferRange.");
        }

        const std::uint64_t tileColumns =
            (static_cast<std::uint64_t>(imageWidth) + options_.tileWidth - 1U) /
            options_.tileWidth;
        const std::uint64_t tileRows =
            (static_cast<std::uint64_t>(imageHeight) + options_.tileHeight - 1U) /
            options_.tileHeight;
        const std::uint64_t batchesPerTile =
            (static_cast<std::uint64_t>(totalSamples) + options_.samplesPerBatch - 1U) /
            options_.samplesPerBatch;
        if (tileColumns > std::numeric_limits<std::uint64_t>::max() / tileRows ||
            tileColumns * tileRows >
                std::numeric_limits<std::uint64_t>::max() / batchesPerTile) {
            throw std::invalid_argument(
                "The configured path render has too many tile dispatches.");
        }
    }

    [[nodiscard]] PathRenderPushConstants makePushConstants(
        const RenderSettings& settings,
        std::uint32_t imageWidth,
        std::uint32_t imageHeight,
        std::uint32_t tileX,
        std::uint32_t tileY,
        std::uint32_t tileWidth,
        std::uint32_t tileHeight,
        std::uint32_t sampleBase,
        std::uint32_t batchCount,
        std::uint32_t totalSamples) const noexcept {
        PathRenderPushConstants result;
        result.cameraPositionAndPixelScale = {
            settings.position.x,
            settings.position.y,
            settings.position.z,
            settings.pixelScale,
        };
        const float rayMinimum =
            std::nextafter(kRayEpsilon, std::numeric_limits<float>::infinity());
        result.cameraDirectionAndRayMinimum = {
            settings.direction.x,
            settings.direction.y,
            settings.direction.z,
            rayMinimum,
        };
        result.cameraUpAndDirectProbability = {
            settings.up.x,
            settings.up.y,
            settings.up.z,
            settings.directLightProbability,
        };
        result.cameraRightAndSampleWeight = {
            settings.right.x,
            settings.right.y,
            settings.right.z,
            1.0F / static_cast<float>(totalSamples),
        };
        result.imageExtentAndTileOrigin = {
            imageWidth,
            imageHeight,
            tileX,
            tileY,
        };
        result.tileExtentAndSampleRange = {
            tileWidth,
            tileHeight,
            sampleBase,
            batchCount,
        };
        result.seedLightCountAndEnvironmentExtent = {
            settings.seed,
            areaLightCount_,
            environmentWidth_,
            environmentHeight_,
        };
        result.environmentAndSceneFlags = {
            environmentFlags_,
            totalSamples,
            sceneVersion_,
            runtime_.primitiveRemapRequired() ? 1U : 0U,
        };
        return result;
    }

    void renderTile(
        const std::uint32_t queueIndex,
        QueueState& queueState,
        const RenderSettings& settings,
        std::uint32_t imageWidth,
        std::uint32_t imageHeight,
        std::uint32_t tileX,
        std::uint32_t tileY,
        std::uint32_t tileWidth,
        std::uint32_t tileHeight,
        std::uint32_t totalSamples,
        Film& film,
        std::uint64_t& directLightSamples,
        double& gpuDispatchMilliseconds,
        std::uint64_t& dispatchCount) {
        const std::size_t tilePixelCount =
            static_cast<std::size_t>(tileWidth) * static_cast<std::size_t>(tileHeight);
        const VkDeviceSize outputBytes = detail::checkedVulkanByteSize(
            tilePixelCount,
            sizeof(PathOutputPixel),
            "Vulkan path-render tile output");
        if (tilePixelCount > queueState.outputCapacity ||
            !queueState.outputBuffer || !queueState.readbackBuffer) {
            throw std::logic_error(
                "A Vulkan path-render queue has insufficient output capacity.");
        }

        const VkPhysicalDeviceLimits& limits = runtime_.physicalProperties().limits;
        const std::uint32_t groupCountX = checkedGroupCount(
            tileWidth,
            kLocalSizeX,
            limits.maxComputeWorkGroupCount[0],
            "X");
        const std::uint32_t groupCountY = checkedGroupCount(
            tileHeight,
            kLocalSizeY,
            limits.maxComputeWorkGroupCount[1],
            "Y");

        VkCommandBuffer commandBuffer = runtime_.beginCommands(queueIndex);
        if (timestampQueriesAvailable_) {
            vkCmdResetQueryPool(
                commandBuffer,
                queueState.timestampQueryPool.get(),
                0U,
                kTimestampQueryCount);
        }
        vkCmdBindPipeline(
            commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_.get());
        vkCmdBindDescriptorSets(
            commandBuffer,
            VK_PIPELINE_BIND_POINT_COMPUTE,
            pipelineLayout_.get(),
            0U,
            1U,
            &queueState.descriptorSet,
            0U,
            nullptr);
        if (timestampQueriesAvailable_) {
            vkCmdWriteTimestamp(
                commandBuffer,
                VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT,
                queueState.timestampQueryPool.get(),
                0U);
        }

        for (std::uint32_t sampleBase = 0U; sampleBase < totalSamples;) {
            const std::uint32_t batchCount =
                std::min(options_.samplesPerBatch, totalSamples - sampleBase);
            if (sampleBase > 0U) {
                const VkMemoryBarrier previousBatchBarrier{
                    VK_STRUCTURE_TYPE_MEMORY_BARRIER,
                    nullptr,
                    VK_ACCESS_SHADER_WRITE_BIT,
                    VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                };
                vkCmdPipelineBarrier(
                    commandBuffer,
                    VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                    VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                    0U,
                    1U,
                    &previousBatchBarrier,
                    0U,
                    nullptr,
                    0U,
                    nullptr);
            }
            const PathRenderPushConstants pushConstants = makePushConstants(
                settings,
                imageWidth,
                imageHeight,
                tileX,
                tileY,
                tileWidth,
                tileHeight,
                sampleBase,
                batchCount,
                totalSamples);
            vkCmdPushConstants(
                commandBuffer,
                pipelineLayout_.get(),
                VK_SHADER_STAGE_COMPUTE_BIT,
                0U,
                sizeof(pushConstants),
                &pushConstants);
            vkCmdDispatch(commandBuffer, groupCountX, groupCountY, 1U);
            ++dispatchCount;
            sampleBase += batchCount;
        }

        if (timestampQueriesAvailable_) {
            vkCmdWriteTimestamp(
                commandBuffer,
                VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT,
                queueState.timestampQueryPool.get(),
                1U);
        }
        const VkBufferMemoryBarrier computeToCopyBarrier{
            VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER,
            nullptr,
            VK_ACCESS_SHADER_WRITE_BIT,
            VK_ACCESS_TRANSFER_READ_BIT,
            VK_QUEUE_FAMILY_IGNORED,
            VK_QUEUE_FAMILY_IGNORED,
            queueState.outputBuffer->get(),
            0U,
            outputBytes,
        };
        vkCmdPipelineBarrier(
            commandBuffer,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            0U,
            0U,
            nullptr,
            1U,
            &computeToCopyBarrier,
            0U,
            nullptr);
        const VkBufferCopy copyRegion{0U, 0U, outputBytes};
        vkCmdCopyBuffer(
            commandBuffer,
            queueState.outputBuffer->get(),
            queueState.readbackBuffer->get(),
            1U,
            &copyRegion);
        const VkBufferMemoryBarrier copyToHostBarrier{
            VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER,
            nullptr,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_ACCESS_HOST_READ_BIT,
            VK_QUEUE_FAMILY_IGNORED,
            VK_QUEUE_FAMILY_IGNORED,
            queueState.readbackBuffer->get(),
            0U,
            outputBytes,
        };
        vkCmdPipelineBarrier(
            commandBuffer,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VK_PIPELINE_STAGE_HOST_BIT,
            0U,
            0U,
            nullptr,
            1U,
            &copyToHostBarrier,
            0U,
            nullptr);
        runtime_.submitAndWait("Vulkan path-render tile", queueIndex);
        if (timestampQueriesAvailable_) {
            gpuDispatchMilliseconds += readGpuTimestampMilliseconds(queueState);
        }
        queueState.readbackBuffer->read(queueState.hostPixels.data(), outputBytes);
        copyTileToFilm(
            queueState.hostPixels,
            tileX,
            tileY,
            tileWidth,
            tileHeight,
            imageWidth,
            settings.exposure,
            film,
            directLightSamples);
    }

    void copyTileToFilm(
        const std::vector<PathOutputPixel>& hostPixels,
        std::uint32_t tileX,
        std::uint32_t tileY,
        std::uint32_t tileWidth,
        std::uint32_t tileHeight,
        std::uint32_t imageWidth,
        float exposure,
        Film& film,
        std::uint64_t& directLightSamples) const {
        for (std::uint32_t localY = 0U; localY < tileHeight; ++localY) {
            for (std::uint32_t localX = 0U; localX < tileWidth; ++localX) {
                const std::size_t localIndex =
                    static_cast<std::size_t>(localY) * tileWidth + localX;
                const std::size_t imageIndex =
                    static_cast<std::size_t>(tileY + localY) * imageWidth +
                    (tileX + localX);
                const PathOutputPixel& source = hostPixels[localIndex];
                const std::uint32_t entering =
                    source.materialEnteringHitAndDirectSampleCount[1];
                const std::uint32_t hit =
                    source.materialEnteringHitAndDirectSampleCount[2];
                if (entering > 1U || hit > 1U) {
                    throw std::runtime_error(
                        "The path-render shader returned an invalid Boolean field.");
                }

                HitInfo& geometry = film.gBuffer[imageIndex];
                geometry.emission = xyz(source.emissionAndEta);
                if (hit != 0U) {
                    const std::uint32_t materialId =
                        source.materialEnteringHitAndDirectSampleCount[0];
                    if (materialId >= materialCount_) {
                        throw std::runtime_error(
                            "The path-render shader returned an invalid material ID.");
                    }
                    geometry.shapeNormal = xyz(source.shapeNormalAndRoughness);
                    geometry.surfaceNormal = xyz(source.surfaceNormalAndSpecular);
                    geometry.position = xyz(source.positionAndMetallic);
                    geometry.baseColor = xyz(source.baseColorAndOpacity);
                    geometry.roughness = source.shapeNormalAndRoughness[3];
                    geometry.specular = source.surfaceNormalAndSpecular[3];
                    geometry.metallic = source.positionAndMetallic[3];
                    geometry.opacity = source.baseColorAndOpacity[3];
                    geometry.eta = source.emissionAndEta[3];
                    geometry.materialId = static_cast<int>(materialId);
                    geometry.entering = entering != 0U;
                }

                assignRadiance(
                    source.directDiffuse,
                    exposure,
                    film.directDiffuseRadiance[imageIndex]);
                assignRadiance(
                    source.directSpecular,
                    exposure,
                    film.directSpecularRadiance[imageIndex]);
                assignRadiance(
                    source.indirectDiffuse,
                    exposure,
                    film.indirectDiffuseRadiance[imageIndex]);
                assignRadiance(
                    source.indirectSpecular,
                    exposure,
                    film.indirectSpecularRadiance[imageIndex]);

                const std::uint64_t sampleCount =
                    source.materialEnteringHitAndDirectSampleCount[3];
                if (sampleCount >
                    std::numeric_limits<std::uint64_t>::max() - directLightSamples) {
                    throw std::overflow_error(
                        "The GPU direct-light sample count overflowed.");
                }
                directLightSamples += sampleCount;
            }
        }
    }

    void createPipeline(const std::vector<std::uint32_t>& shaderCode) {
        const VkDevice device = runtime_.device();
        std::array<VkDescriptorSetLayoutBinding, kStorageBufferBindingCount + 1U>
            bindings{};
        bindings[0] = VkDescriptorSetLayoutBinding{
            0U,
            VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR,
            1U,
            VK_SHADER_STAGE_COMPUTE_BIT,
            nullptr,
        };
        for (std::uint32_t binding = 1U; binding < bindings.size(); ++binding) {
            bindings[binding] = VkDescriptorSetLayoutBinding{
                binding,
                VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
                1U,
                VK_SHADER_STAGE_COMPUTE_BIT,
                nullptr,
            };
        }
        const VkDescriptorSetLayoutCreateInfo descriptorLayoutInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,
            nullptr,
            0U,
            static_cast<std::uint32_t>(bindings.size()),
            bindings.data(),
        };
        VkDescriptorSetLayout descriptorLayout = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateDescriptorSetLayout(
                device, &descriptorLayoutInfo, nullptr, &descriptorLayout),
            "vkCreateDescriptorSetLayout(path render)");
        descriptorLayout_.reset(
            descriptorLayout,
            [device](VkDescriptorSetLayout value) {
                vkDestroyDescriptorSetLayout(device, value, nullptr);
            });

        const VkPushConstantRange pushConstantRange{
            VK_SHADER_STAGE_COMPUTE_BIT,
            0U,
            sizeof(PathRenderPushConstants),
        };
        const VkDescriptorSetLayout descriptorLayoutHandle = descriptorLayout_.get();
        const VkPipelineLayoutCreateInfo pipelineLayoutInfo{
            VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,
            nullptr,
            0U,
            1U,
            &descriptorLayoutHandle,
            1U,
            &pushConstantRange,
        };
        VkPipelineLayout pipelineLayout = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreatePipelineLayout(
                device, &pipelineLayoutInfo, nullptr, &pipelineLayout),
            "vkCreatePipelineLayout(path render)");
        pipelineLayout_.reset(
            pipelineLayout,
            [device](VkPipelineLayout value) {
                vkDestroyPipelineLayout(device, value, nullptr);
            });

        const VkShaderModuleCreateInfo shaderModuleInfo{
            VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,
            nullptr,
            0U,
            shaderCode.size() * sizeof(std::uint32_t),
            shaderCode.data(),
        };
        VkShaderModule shaderModule = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateShaderModule(device, &shaderModuleInfo, nullptr, &shaderModule),
            "vkCreateShaderModule(path render)");
        shaderModule_.reset(
            shaderModule,
            [device](VkShaderModule value) {
                vkDestroyShaderModule(device, value, nullptr);
            });

        const VkPipelineShaderStageCreateInfo shaderStageInfo{
            VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
            nullptr,
            0U,
            VK_SHADER_STAGE_COMPUTE_BIT,
            shaderModule_.get(),
            "main",
            nullptr,
        };
        const VkComputePipelineCreateInfo pipelineInfo{
            VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO,
            nullptr,
            0U,
            shaderStageInfo,
            pipelineLayout_.get(),
            VK_NULL_HANDLE,
            0,
        };
        VkPipeline pipeline = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateComputePipelines(
                device, VK_NULL_HANDLE, 1U, &pipelineInfo, nullptr, &pipeline),
            "vkCreateComputePipelines(path render)");
        pipeline_.reset(
            pipeline,
            [device](VkPipeline value) {
                vkDestroyPipeline(device, value, nullptr);
            });

        const std::uint32_t queueCount = runtime_.computeQueueCount();
        if (queueCount >
            std::numeric_limits<std::uint32_t>::max() /
                kStorageBufferBindingCount) {
            throw std::overflow_error(
                "The Vulkan path-render descriptor count overflowed.");
        }
        const std::array<VkDescriptorPoolSize, 2> poolSizes{{
            {VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR, queueCount},
            {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
             queueCount * kStorageBufferBindingCount},
        }};
        const VkDescriptorPoolCreateInfo descriptorPoolInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
            nullptr,
            0U,
            queueCount,
            static_cast<std::uint32_t>(poolSizes.size()),
            poolSizes.data(),
        };
        VkDescriptorPool descriptorPool = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateDescriptorPool(device, &descriptorPoolInfo, nullptr, &descriptorPool),
            "vkCreateDescriptorPool(path render)");
        descriptorPool_.reset(
            descriptorPool,
            [device](VkDescriptorPool value) {
                vkDestroyDescriptorPool(device, value, nullptr);
            });
        const std::vector<VkDescriptorSetLayout> descriptorLayouts(
            queueCount, descriptorLayoutHandle);
        std::vector<VkDescriptorSet> descriptorSets(queueCount, VK_NULL_HANDLE);
        const VkDescriptorSetAllocateInfo allocateInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
            nullptr,
            descriptorPool_.get(),
            queueCount,
            descriptorLayouts.data(),
        };
        detail::requireVulkanSuccess(
            vkAllocateDescriptorSets(device, &allocateInfo, descriptorSets.data()),
            "vkAllocateDescriptorSets(path render)");
        for (std::uint32_t queueIndex = 0U;
             queueIndex < queueCount;
             ++queueIndex) {
            queueStates_[queueIndex].descriptorSet = descriptorSets[queueIndex];
        }
    }

    void createTimestampQueries() {
        const std::uint32_t validBits = runtime_.timestampValidBits();
        const float timestampPeriod = runtime_.physicalProperties().limits.timestampPeriod;
        timestampQueriesAvailable_ =
            validBits != 0U && std::isfinite(timestampPeriod) && timestampPeriod > 0.0F;
        if (!timestampQueriesAvailable_) {
            return;
        }

        const VkQueryPoolCreateInfo queryPoolInfo{
            VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO,
            nullptr,
            0U,
            VK_QUERY_TYPE_TIMESTAMP,
            kTimestampQueryCount,
            0U,
        };
        const VkDevice device = runtime_.device();
        for (QueueState& queueState : queueStates_) {
            VkQueryPool queryPool = VK_NULL_HANDLE;
            detail::requireVulkanSuccess(
                vkCreateQueryPool(device, &queryPoolInfo, nullptr, &queryPool),
                "vkCreateQueryPool(path timestamp)");
            queueState.timestampQueryPool.reset(
                queryPool,
                [device](VkQueryPool value) {
                    vkDestroyQueryPool(device, value, nullptr);
                });
        }
    }

    void ensureOutputCapacity(
        QueueState& queueState,
        const std::size_t pixelCount,
        const VkDeviceSize outputBytes) {
        if (pixelCount <= queueState.outputCapacity) {
            return;
        }

        constexpr VkMemoryPropertyFlags kHostReadbackMemory =
            VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
        const VkDevice device = runtime_.device();
        auto newOutputBuffer = std::make_unique<detail::VulkanBuffer>(
            device,
            runtime_.memoryProperties(),
            outputBytes,
            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
            false);
        auto newReadbackBuffer = std::make_unique<detail::VulkanBuffer>(
            device,
            runtime_.memoryProperties(),
            outputBytes,
            VK_BUFFER_USAGE_TRANSFER_DST_BIT,
            kHostReadbackMemory,
            false);

        updateDescriptorSet(queueState, *newOutputBuffer);
        queueState.outputBuffer = std::move(newOutputBuffer);
        queueState.readbackBuffer = std::move(newReadbackBuffer);
        queueState.outputCapacity = pixelCount;
    }

    void updateDescriptorSet(
        const QueueState& queueState,
        const detail::VulkanBuffer& outputBuffer) {
        const VkAccelerationStructureKHR topLevel =
            runtime_.topLevelAccelerationStructure();
        const VkWriteDescriptorSetAccelerationStructureKHR accelerationWrite{
            VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET_ACCELERATION_STRUCTURE_KHR,
            nullptr,
            1U,
            &topLevel,
        };
        const std::array<VkDescriptorBufferInfo, kStorageBufferBindingCount>
            bufferInfos{{
                {runtime_.vertexBuffer().get(), 0U, runtime_.vertexBuffer().size()},
                {runtime_.indexBuffer().get(), 0U, runtime_.indexBuffer().size()},
                {runtime_.triangleMaterialIdBuffer().get(),
                 0U,
                 runtime_.triangleMaterialIdBuffer().size()},
                {runtime_.materialBuffer().get(), 0U, runtime_.materialBuffer().size()},
                {outputBuffer.get(), 0U, outputBuffer.size()},
                {runtime_.textureDescriptorBuffer().get(),
                 0U,
                 runtime_.textureDescriptorBuffer().size()},
                {runtime_.textureMipBuffer().get(),
                 0U,
                 runtime_.textureMipBuffer().size()},
                {runtime_.textureTexelBuffer().get(),
                 0U,
                 runtime_.textureTexelBuffer().size()},
                {runtime_.areaLightBuffer().get(),
                 0U,
                 runtime_.areaLightBuffer().size()},
                {runtime_.areaLightTriangleBuffer().get(),
                 0U,
                 runtime_.areaLightTriangleBuffer().size()},
                {runtime_.environmentRowBuffer().get(),
                 0U,
                 runtime_.environmentRowBuffer().size()},
                {runtime_.environmentTexelBuffer().get(),
                 0U,
                 runtime_.environmentTexelBuffer().size()},
                {runtime_.primitiveRemapBuffer().get(),
                 0U,
                 runtime_.primitiveRemapBuffer().size()},
            }};

        std::array<VkWriteDescriptorSet, kStorageBufferBindingCount + 1U> writes{};
        writes[0].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[0].pNext = &accelerationWrite;
        writes[0].dstSet = queueState.descriptorSet;
        writes[0].dstBinding = 0U;
        writes[0].descriptorCount = 1U;
        writes[0].descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
        for (std::uint32_t binding = 1U; binding < writes.size(); ++binding) {
            VkWriteDescriptorSet& write = writes[binding];
            write.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            write.dstSet = queueState.descriptorSet;
            write.dstBinding = binding;
            write.descriptorCount = 1U;
            write.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            write.pBufferInfo = &bufferInfos[binding - 1U];
        }
        vkUpdateDescriptorSets(
            runtime_.device(),
            static_cast<std::uint32_t>(writes.size()),
            writes.data(),
            0U,
            nullptr);
    }

    [[nodiscard]] double readGpuTimestampMilliseconds(
        const QueueState& queueState) const {
        std::array<std::uint64_t, kTimestampQueryCount> timestamps{};
        detail::requireVulkanSuccess(
            vkGetQueryPoolResults(
                runtime_.device(),
                queueState.timestampQueryPool.get(),
                0U,
                kTimestampQueryCount,
                sizeof(timestamps),
                timestamps.data(),
                sizeof(std::uint64_t),
                VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT),
            "vkGetQueryPoolResults(path timestamp)");
        const std::uint64_t ticks = timestampDelta(
            timestamps[0], timestamps[1], runtime_.timestampValidBits());
        const double nanoseconds =
            static_cast<double>(ticks) *
            static_cast<double>(runtime_.physicalProperties().limits.timestampPeriod);
        return nanoseconds / 1.0e6;
    }

    VulkanPathRenderOptions options_;
    // Declared before all renderer-owned Vulkan objects so it is destroyed after them.
    detail::VulkanRuntime runtime_;
    std::uint32_t materialCount_{0U};
    std::uint32_t areaLightCount_{0U};
    std::uint32_t environmentWidth_{0U};
    std::uint32_t environmentHeight_{0U};
    std::uint32_t environmentFlags_{0U};
    std::uint32_t sceneVersion_{0U};
    bool timestampQueriesAvailable_{false};
    detail::UniqueVulkanHandle<VkDescriptorSetLayout> descriptorLayout_;
    detail::UniqueVulkanHandle<VkPipelineLayout> pipelineLayout_;
    detail::UniqueVulkanHandle<VkShaderModule> shaderModule_;
    detail::UniqueVulkanHandle<VkPipeline> pipeline_;
    detail::UniqueVulkanHandle<VkDescriptorPool> descriptorPool_;
    std::vector<QueueState> queueStates_;
};

VulkanPathRenderer::VulkanPathRenderer(
    const PackedSceneData& scene,
    const std::filesystem::path& spirvPath,
    VulkanPathRenderOptions options)
    : implementation_(
          std::make_unique<Implementation>(scene, spirvPath, std::move(options))) {}

VulkanPathRenderer::~VulkanPathRenderer() = default;

const std::string& VulkanPathRenderer::deviceName() const noexcept {
    return implementation_->deviceName();
}

std::uint32_t VulkanPathRenderer::computeQueueCount() const noexcept {
    return implementation_->computeQueueCount();
}

const VulkanRayQuerySetupTimings& VulkanPathRenderer::setupTimings() const noexcept {
    return implementation_->setupTimings();
}

VulkanValidationReport VulkanPathRenderer::validationReport() const {
    return implementation_->validationReport();
}

VulkanPathRenderResult VulkanPathRenderer::render(const RenderSettings& settings) {
    return implementation_->render(settings);
}

}  // namespace raym0nade::gpu
