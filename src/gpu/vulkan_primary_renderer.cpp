#include "raym0nade/gpu/vulkan_primary_renderer.hpp"

#include <vulkan/vulkan.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
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

struct alignas(16) PrimaryRenderPushConstants {
    std::array<float, 4> cameraPositionAndPixelScale{};
    std::array<float, 4> cameraDirectionAndRayMinimum{};
    std::array<float, 4> cameraUpAndUnused{};
    std::array<float, 4> cameraRightAndUnused{};
    std::array<std::uint32_t, 4> extentAovAndReserved{};
    std::array<float, 4> directionToLightAndUnused{};
    std::array<float, 4> incidentRadianceAndUnused{};
};

struct alignas(16) PrimaryOutputPixel {
    std::array<float, 4> color{};
};

static_assert(sizeof(PrimaryRenderPushConstants) == 112U);
static_assert(alignof(PrimaryRenderPushConstants) == 16U);
static_assert(offsetof(PrimaryRenderPushConstants, cameraPositionAndPixelScale) == 0U);
static_assert(offsetof(PrimaryRenderPushConstants, cameraDirectionAndRayMinimum) == 16U);
static_assert(offsetof(PrimaryRenderPushConstants, cameraUpAndUnused) == 32U);
static_assert(offsetof(PrimaryRenderPushConstants, cameraRightAndUnused) == 48U);
static_assert(offsetof(PrimaryRenderPushConstants, extentAovAndReserved) == 64U);
static_assert(offsetof(PrimaryRenderPushConstants, directionToLightAndUnused) == 80U);
static_assert(offsetof(PrimaryRenderPushConstants, incidentRadianceAndUnused) == 96U);
static_assert(std::is_standard_layout_v<PrimaryRenderPushConstants>);
static_assert(std::is_trivially_copyable_v<PrimaryRenderPushConstants>);
static_assert(sizeof(PrimaryOutputPixel) == 16U);
static_assert(alignof(PrimaryOutputPixel) == 16U);
static_assert(std::is_standard_layout_v<PrimaryOutputPixel>);
static_assert(std::is_trivially_copyable_v<PrimaryOutputPixel>);
static_assert(static_cast<std::uint32_t>(PrimaryAov::BaseColor) == 0U);
static_assert(static_cast<std::uint32_t>(PrimaryAov::ShapeNormal) == 1U);
static_assert(static_cast<std::uint32_t>(PrimaryAov::DirectDiffuse) == 2U);

[[nodiscard]] bool referencedMaterialsHaveFlag(
    const PackedSceneData& scene, std::uint32_t flag) noexcept {
    return std::any_of(
        scene.triangleMaterialIds.begin(),
        scene.triangleMaterialIds.end(),
        [&scene, flag](std::uint32_t materialId) {
            return (scene.materials[materialId].flagsAndReserved[0] & flag) != 0U;
        });
}

[[nodiscard]] bool referencedMaterialsAreNotOpaque(const PackedSceneData& scene) noexcept {
    return std::any_of(
        scene.triangleMaterialIds.begin(),
        scene.triangleMaterialIds.end(),
        [&scene](std::uint32_t materialId) {
            return scene.materials[materialId].diffuseAndOpacity[3] <=
                   1.0F - kRayEpsilon;
        });
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
            std::string{"The primary render exceeds the Vulkan "} + axis +
            " dispatch limit.");
    }
    return static_cast<std::uint32_t>(groups);
}

[[nodiscard]] PrimaryRenderPushConstants makePushConstants(
    const PrimaryRenderRequest& request) noexcept {
    PrimaryRenderPushConstants result;
    const vec3 directionToLight =
        safeNormalize(request.directionalLight.directionToLight);
    result.cameraPositionAndPixelScale = {
        request.camera.position.x,
        request.camera.position.y,
        request.camera.position.z,
        request.camera.pixelScale,
    };
    // Vulkan ray queries include the lower interval endpoint, while the CPU triangle test uses a
    // strict distance > kRayEpsilon comparison. Advancing by one float makes the intervals match.
    const float rayMinimum =
        std::nextafter(kRayEpsilon, std::numeric_limits<float>::infinity());
    result.cameraDirectionAndRayMinimum = {
        request.camera.direction.x,
        request.camera.direction.y,
        request.camera.direction.z,
        rayMinimum,
    };
    result.cameraUpAndUnused = {
        request.camera.up.x,
        request.camera.up.y,
        request.camera.up.z,
        0.0F,
    };
    result.cameraRightAndUnused = {
        request.camera.right.x,
        request.camera.right.y,
        request.camera.right.z,
        0.0F,
    };
    result.extentAovAndReserved = {
        request.extent.width,
        request.extent.height,
        static_cast<std::uint32_t>(request.aov),
        0U,
    };
    result.directionToLightAndUnused = {
        directionToLight.x,
        directionToLight.y,
        directionToLight.z,
        0.0F,
    };
    result.incidentRadianceAndUnused = {
        request.directionalLight.incidentRadiance.x,
        request.directionalLight.incidentRadiance.y,
        request.directionalLight.incidentRadiance.z,
        0.0F,
    };
    return result;
}

[[nodiscard]] std::uint64_t timestampDelta(
    std::uint64_t begin, std::uint64_t end, std::uint32_t validBits) noexcept {
    if (validBits >= 64U) {
        return end - begin;
    }
    const std::uint64_t mask = (std::uint64_t{1U} << validBits) - 1U;
    return (end - begin) & mask;
}

}  // namespace

class VulkanPrimaryRenderer::Implementation {
public:
    Implementation(
        const PackedSceneData& scene,
        const std::filesystem::path& spirvPath,
        VulkanRayQueryOptions options)
        : runtime_(scene, options),
          hasReferencedDiffuseTexture_(
              referencedMaterialsHaveFlag(scene, kPackedMaterialHasDiffuseTexture)),
          hasReferencedTransparentMaterial_(referencedMaterialsAreNotOpaque(scene)) {
        validateDeviceLimits();
        createPipeline(detail::readSpirvFile(spirvPath));
        createTimestampQueries();
    }

    [[nodiscard]] const std::string& deviceName() const noexcept {
        return runtime_.deviceName();
    }

    [[nodiscard]] const VulkanRayQuerySetupTimings& setupTimings() const noexcept {
        return runtime_.setupTimings();
    }

    [[nodiscard]] VulkanValidationReport validationReport() const {
        return runtime_.validationReport();
    }

    [[nodiscard]] VulkanPrimaryRenderResult render(const PrimaryRenderRequest& request) {
        request.validate();
        validateAovFeatures(request.aov);

        const std::size_t pixelCount = request.extent.pixelCount();
        const VkDeviceSize outputBytes = detail::checkedVulkanByteSize(
            pixelCount, sizeof(PrimaryOutputPixel), "Vulkan primary render output");
        const VkPhysicalDeviceLimits& limits = runtime_.physicalProperties().limits;
        if (outputBytes > limits.maxStorageBufferRange) {
            throw std::invalid_argument(
                "The primary render output exceeds maxStorageBufferRange.");
        }
        const std::uint32_t groupCountX = checkedGroupCount(
            request.extent.width,
            kLocalSizeX,
            limits.maxComputeWorkGroupCount[0],
            "X");
        const std::uint32_t groupCountY = checkedGroupCount(
            request.extent.height,
            kLocalSizeY,
            limits.maxComputeWorkGroupCount[1],
            "Y");

        std::lock_guard<std::mutex> lock{runtime_.operationMutex()};
        ensureOutputCapacity(pixelCount, outputBytes);
        hostPixels_.resize(pixelCount);

        VulkanPrimaryRenderResult result;
        result.image.extent = request.extent;
        result.image.pixels.resize(pixelCount);
        result.timings.gpuTimestampAvailable = timestampQueriesAvailable_;

        const detail::VulkanClock::time_point begin = detail::VulkanClock::now();
        VkCommandBuffer commandBuffer = runtime_.beginCommands();
        if (timestampQueriesAvailable_) {
            vkCmdResetQueryPool(commandBuffer, timestampQueryPool_.get(), 0U, kTimestampQueryCount);
        }
        vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_.get());
        vkCmdBindDescriptorSets(
            commandBuffer,
            VK_PIPELINE_BIND_POINT_COMPUTE,
            pipelineLayout_.get(),
            0U,
            1U,
            &descriptorSet_,
            0U,
            nullptr);
        const PrimaryRenderPushConstants pushConstants = makePushConstants(request);
        vkCmdPushConstants(
            commandBuffer,
            pipelineLayout_.get(),
            VK_SHADER_STAGE_COMPUTE_BIT,
            0U,
            sizeof(pushConstants),
            &pushConstants);
        if (timestampQueriesAvailable_) {
            vkCmdWriteTimestamp(
                commandBuffer,
                VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT,
                timestampQueryPool_.get(),
                0U);
        }
        vkCmdDispatch(commandBuffer, groupCountX, groupCountY, 1U);
        if (timestampQueriesAvailable_) {
            vkCmdWriteTimestamp(
                commandBuffer,
                VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT,
                timestampQueryPool_.get(),
                1U);
        }

        const VkMemoryBarrier computeToCopyBarrier{
            VK_STRUCTURE_TYPE_MEMORY_BARRIER,
            nullptr,
            VK_ACCESS_SHADER_WRITE_BIT,
            VK_ACCESS_TRANSFER_READ_BIT,
        };
        vkCmdPipelineBarrier(
            commandBuffer,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            0U,
            1U,
            &computeToCopyBarrier,
            0U,
            nullptr,
            0U,
            nullptr);
        const VkBufferCopy copyRegion{0U, 0U, outputBytes};
        vkCmdCopyBuffer(
            commandBuffer,
            outputBuffer_->get(),
            readbackBuffer_->get(),
            1U,
            &copyRegion);
        const VkMemoryBarrier copyToHostBarrier{
            VK_STRUCTURE_TYPE_MEMORY_BARRIER,
            nullptr,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_ACCESS_HOST_READ_BIT,
        };
        vkCmdPipelineBarrier(
            commandBuffer,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VK_PIPELINE_STAGE_HOST_BIT,
            0U,
            1U,
            &copyToHostBarrier,
            0U,
            nullptr,
            0U,
            nullptr);
        runtime_.submitAndWait("Vulkan primary-AOV render");

        readbackBuffer_->read(hostPixels_.data(), outputBytes);
        result.timings.dispatchAndReadbackMilliseconds =
            detail::elapsedMilliseconds(begin);
        readGpuTimestamps(result.timings);

        for (std::size_t index = 0U; index < pixelCount; ++index) {
            const std::array<float, 4>& color = hostPixels_[index].color;
            result.image.pixels[index] = vec3{color[0], color[1], color[2]};
        }
        result.image.validate();
        return result;
    }

private:
    void validateDeviceLimits() const {
        const VkPhysicalDeviceLimits& limits = runtime_.physicalProperties().limits;
        if (limits.maxComputeWorkGroupSize[0] < kLocalSizeX ||
            limits.maxComputeWorkGroupSize[1] < kLocalSizeY ||
            limits.maxComputeWorkGroupInvocations < kLocalSizeX * kLocalSizeY) {
            throw std::runtime_error(
                "The Vulkan device cannot execute the primary renderer's 8 x 8 workgroup.");
        }
        if (limits.maxPushConstantsSize < sizeof(PrimaryRenderPushConstants)) {
            throw std::runtime_error(
                "The Vulkan device does not provide enough push-constant storage.");
        }
        constexpr std::uint32_t kStorageBufferBindingCount = 5U;
        if (limits.maxPerStageDescriptorStorageBuffers < kStorageBufferBindingCount ||
            limits.maxDescriptorSetStorageBuffers < kStorageBufferBindingCount) {
            throw std::runtime_error(
                "The Vulkan device does not provide enough storage-buffer descriptors.");
        }
        if (runtime_.vertexBuffer().size() > limits.maxStorageBufferRange ||
            runtime_.indexBuffer().size() > limits.maxStorageBufferRange ||
            runtime_.triangleMaterialIdBuffer().size() > limits.maxStorageBufferRange ||
            runtime_.materialBuffer().size() > limits.maxStorageBufferRange) {
            throw std::runtime_error(
                "A packed-scene buffer exceeds the Vulkan storage-buffer range limit.");
        }
    }

    void validateAovFeatures(PrimaryAov aov) const {
        if (aov != PrimaryAov::BaseColor && aov != PrimaryAov::DirectDiffuse) {
            return;
        }
        const char* aovName =
            aov == PrimaryAov::BaseColor ? "BaseColor" : "DirectDiffuse";
        if (hasReferencedDiffuseTexture_) {
            throw std::invalid_argument(
                std::string{"Vulkan "} + aovName +
                " rendering does not yet support diffuse textures.");
        }
        if (hasReferencedTransparentMaterial_) {
            throw std::invalid_argument(
                std::string{"Vulkan "} + aovName +
                " rendering currently requires fully opaque materials.");
        }
    }

    void createPipeline(const std::vector<std::uint32_t>& shaderCode) {
        const VkDevice device = runtime_.device();
        std::array<VkDescriptorSetLayoutBinding, 6> bindings{};
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
            "vkCreateDescriptorSetLayout");
        descriptorLayout_.reset(
            descriptorLayout,
            [device](VkDescriptorSetLayout value) {
                vkDestroyDescriptorSetLayout(device, value, nullptr);
            });

        const VkPushConstantRange pushConstantRange{
            VK_SHADER_STAGE_COMPUTE_BIT,
            0U,
            sizeof(PrimaryRenderPushConstants),
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
            vkCreatePipelineLayout(device, &pipelineLayoutInfo, nullptr, &pipelineLayout),
            "vkCreatePipelineLayout");
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
            "vkCreateShaderModule");
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
            "vkCreateComputePipelines");
        pipeline_.reset(
            pipeline,
            [device](VkPipeline value) { vkDestroyPipeline(device, value, nullptr); });

        const std::array<VkDescriptorPoolSize, 2> poolSizes{{
            {VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR, 1U},
            {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 5U},
        }};
        const VkDescriptorPoolCreateInfo descriptorPoolInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
            nullptr,
            0U,
            1U,
            static_cast<std::uint32_t>(poolSizes.size()),
            poolSizes.data(),
        };
        VkDescriptorPool descriptorPool = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateDescriptorPool(device, &descriptorPoolInfo, nullptr, &descriptorPool),
            "vkCreateDescriptorPool");
        descriptorPool_.reset(
            descriptorPool,
            [device](VkDescriptorPool value) {
                vkDestroyDescriptorPool(device, value, nullptr);
            });
        const VkDescriptorSetAllocateInfo allocateInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
            nullptr,
            descriptorPool_.get(),
            1U,
            &descriptorLayoutHandle,
        };
        detail::requireVulkanSuccess(
            vkAllocateDescriptorSets(device, &allocateInfo, &descriptorSet_),
            "vkAllocateDescriptorSets");
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
        VkQueryPool queryPool = VK_NULL_HANDLE;
        const VkDevice device = runtime_.device();
        detail::requireVulkanSuccess(
            vkCreateQueryPool(device, &queryPoolInfo, nullptr, &queryPool),
            "vkCreateQueryPool(timestamp)");
        timestampQueryPool_.reset(
            queryPool,
            [device](VkQueryPool value) { vkDestroyQueryPool(device, value, nullptr); });
    }

    void ensureOutputCapacity(std::size_t pixelCount, VkDeviceSize outputBytes) {
        if (pixelCount <= outputCapacity_) {
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

        updateDescriptorSet(*newOutputBuffer);
        outputBuffer_ = std::move(newOutputBuffer);
        readbackBuffer_ = std::move(newReadbackBuffer);
        outputCapacity_ = pixelCount;
    }

    void updateDescriptorSet(const detail::VulkanBuffer& outputBuffer) {
        const VkAccelerationStructureKHR topLevel =
            runtime_.topLevelAccelerationStructure();
        const VkWriteDescriptorSetAccelerationStructureKHR accelerationWrite{
            VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET_ACCELERATION_STRUCTURE_KHR,
            nullptr,
            1U,
            &topLevel,
        };
        const std::array<VkDescriptorBufferInfo, 5> bufferInfos{{
            {runtime_.vertexBuffer().get(), 0U, runtime_.vertexBuffer().size()},
            {runtime_.indexBuffer().get(), 0U, runtime_.indexBuffer().size()},
            {runtime_.triangleMaterialIdBuffer().get(),
             0U,
             runtime_.triangleMaterialIdBuffer().size()},
            {runtime_.materialBuffer().get(), 0U, runtime_.materialBuffer().size()},
            {outputBuffer.get(), 0U, outputBuffer.size()},
        }};

        std::array<VkWriteDescriptorSet, 6> writes{};
        writes[0].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[0].pNext = &accelerationWrite;
        writes[0].dstSet = descriptorSet_;
        writes[0].dstBinding = 0U;
        writes[0].descriptorCount = 1U;
        writes[0].descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
        for (std::uint32_t binding = 1U; binding < writes.size(); ++binding) {
            VkWriteDescriptorSet& write = writes[binding];
            write.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
            write.dstSet = descriptorSet_;
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

    void readGpuTimestamps(VulkanPrimaryRenderTimings& timings) const {
        if (!timestampQueriesAvailable_) {
            return;
        }
        std::array<std::uint64_t, kTimestampQueryCount> timestamps{};
        detail::requireVulkanSuccess(
            vkGetQueryPoolResults(
                runtime_.device(),
                timestampQueryPool_.get(),
                0U,
                kTimestampQueryCount,
                sizeof(timestamps),
                timestamps.data(),
                sizeof(std::uint64_t),
                VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT),
            "vkGetQueryPoolResults(timestamp)");
        const std::uint64_t ticks = timestampDelta(
            timestamps[0], timestamps[1], runtime_.timestampValidBits());
        const double nanoseconds = static_cast<double>(ticks) *
                                   static_cast<double>(
                                       runtime_.physicalProperties().limits.timestampPeriod);
        timings.gpuDispatchMilliseconds = nanoseconds / 1.0e6;
    }

    // Runtime is declared first so every renderer-owned Vulkan object is destroyed before the
    // device and acceleration structures it depends on.
    detail::VulkanRuntime runtime_;
    bool hasReferencedDiffuseTexture_{false};
    bool hasReferencedTransparentMaterial_{false};
    bool timestampQueriesAvailable_{false};
    std::size_t outputCapacity_{0U};
    std::vector<PrimaryOutputPixel> hostPixels_;
    std::unique_ptr<detail::VulkanBuffer> outputBuffer_;
    std::unique_ptr<detail::VulkanBuffer> readbackBuffer_;
    detail::UniqueVulkanHandle<VkDescriptorSetLayout> descriptorLayout_;
    detail::UniqueVulkanHandle<VkPipelineLayout> pipelineLayout_;
    detail::UniqueVulkanHandle<VkShaderModule> shaderModule_;
    detail::UniqueVulkanHandle<VkPipeline> pipeline_;
    detail::UniqueVulkanHandle<VkDescriptorPool> descriptorPool_;
    VkDescriptorSet descriptorSet_{VK_NULL_HANDLE};
    detail::UniqueVulkanHandle<VkQueryPool> timestampQueryPool_;
};

VulkanPrimaryRenderer::VulkanPrimaryRenderer(
    const PackedSceneData& scene,
    const std::filesystem::path& spirvPath,
    VulkanRayQueryOptions options)
    : implementation_(std::make_unique<Implementation>(scene, spirvPath, options)) {}

VulkanPrimaryRenderer::~VulkanPrimaryRenderer() = default;

const std::string& VulkanPrimaryRenderer::deviceName() const noexcept {
    return implementation_->deviceName();
}

const VulkanRayQuerySetupTimings& VulkanPrimaryRenderer::setupTimings() const noexcept {
    return implementation_->setupTimings();
}

VulkanValidationReport VulkanPrimaryRenderer::validationReport() const {
    return implementation_->validationReport();
}

VulkanPrimaryRenderResult VulkanPrimaryRenderer::render(
    const PrimaryRenderRequest& request) {
    return implementation_->render(request);
}

}  // namespace raym0nade::gpu
