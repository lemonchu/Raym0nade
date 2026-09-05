#include "raym0nade/gpu/vulkan_counter_rng_test.hpp"

#include <vulkan/vulkan.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <mutex>
#include <stdexcept>
#include <vector>

#include "raym0nade/scene_data.hpp"
#include "vulkan_runtime.hpp"

namespace raym0nade::gpu {
namespace {

constexpr std::uint32_t kShaderWorkgroupSize = 4U;
constexpr std::uint32_t kDispatchGroupCount =
    static_cast<std::uint32_t>(
        (kVulkanCounterRngKatAddressCount + kShaderWorkgroupSize - 1U) /
        kShaderWorkgroupSize);

PackedSceneData makeCounterRngKatScene() {
    PackedSceneData scene;
    scene.vertices = {
        PackedVertex{{-1.0F, -1.0F, 0.0F, 0.0F}, {0.0F, 1.0F, 0.0F, 0.0F}},
        PackedVertex{{1.0F, -1.0F, 0.0F, 0.0F}, {0.0F, 1.0F, 1.0F, 0.0F}},
        PackedVertex{{0.0F, 1.0F, 0.0F, 0.0F}, {0.0F, 1.0F, 0.5F, 1.0F}},
    };
    scene.triangleIndices = {0U, 1U, 2U};
    scene.triangleMaterialIds = {0U};
    scene.materials.emplace_back();
    scene.validate();
    return scene;
}

using ObservationArray =
    std::array<VulkanCounterRngObservation, kVulkanCounterRngKatAddressCount>;

ObservationArray makePoisonedObservations(const std::uint32_t salt) noexcept {
    ObservationArray result{};
    for (std::size_t index = 0U; index < result.size(); ++index) {
        const std::uint32_t value = static_cast<std::uint32_t>(index);
        result[index].word = 0xA5A50000U ^ salt ^ value;
        result[index].open01Bits = 0x5A5A0000U ^ salt ^ (value << 8U);
        result[index].blockWord = 0x3C3C0000U ^ salt ^ (value << 12U);
        result[index].blockOpen01Bits = 0xC3C30000U ^ salt ^ (value << 16U);
    }
    return result;
}

class VulkanCounterRngKatRunner {
public:
    VulkanCounterRngKatRunner(
        const std::filesystem::path& spirvPath,
        const VulkanRayQueryOptions options)
        : runtime_(makeCounterRngKatScene(), options),
          outputBuffer_(
              runtime_.device(),
              runtime_.memoryProperties(),
              outputByteSize(),
              VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
              VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                  VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
              false) {
        validateDeviceLimits();
        createPipeline(detail::readSpirvFile(spirvPath));
        createDescriptorSet();
    }

    [[nodiscard]] VulkanCounterRngKatResult run() {
        VulkanCounterRngKatResult result;
        result.deviceName = runtime_.deviceName();
        result.firstDispatch = dispatch(0x13579BDFU);
        result.repeatedDispatch = dispatch(0x2468ACE0U);
        result.validation = runtime_.validationReport();
        return result;
    }

private:
    [[nodiscard]] static VkDeviceSize outputByteSize() {
        return detail::checkedVulkanByteSize(
            kVulkanCounterRngKatAddressCount,
            sizeof(VulkanCounterRngObservation),
            "Vulkan counter RNG KAT output");
    }

    void validateDeviceLimits() const {
        const VkPhysicalDeviceLimits& limits = runtime_.physicalProperties().limits;
        if (outputByteSize() > limits.maxStorageBufferRange) {
            throw std::invalid_argument(
                "The Vulkan counter RNG KAT output exceeds maxStorageBufferRange.");
        }
        if (limits.maxPerStageDescriptorStorageBuffers < 1U ||
            limits.maxDescriptorSetStorageBuffers < 1U) {
            throw std::runtime_error(
                "The Vulkan device cannot bind the counter RNG KAT output buffer.");
        }
        if (limits.maxComputeWorkGroupInvocations < kShaderWorkgroupSize ||
            limits.maxComputeWorkGroupSize[0] < kShaderWorkgroupSize ||
            limits.maxComputeWorkGroupCount[0] < kDispatchGroupCount) {
            throw std::runtime_error(
                "The Vulkan device cannot dispatch the counter RNG KAT workload.");
        }
    }

    void createPipeline(const std::vector<std::uint32_t>& shaderCode) {
        const VkDevice device = runtime_.device();
        const VkDescriptorSetLayoutBinding binding{
            0U,
            VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            1U,
            VK_SHADER_STAGE_COMPUTE_BIT,
            nullptr,
        };
        const VkDescriptorSetLayoutCreateInfo descriptorLayoutInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,
            nullptr,
            0U,
            1U,
            &binding,
        };
        VkDescriptorSetLayout descriptorLayout = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateDescriptorSetLayout(
                device, &descriptorLayoutInfo, nullptr, &descriptorLayout),
            "vkCreateDescriptorSetLayout(counter RNG KAT)");
        descriptorLayout_.reset(
            descriptorLayout,
            [device](const VkDescriptorSetLayout value) {
                vkDestroyDescriptorSetLayout(device, value, nullptr);
            });

        const VkPipelineLayoutCreateInfo pipelineLayoutInfo{
            VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,
            nullptr,
            0U,
            1U,
            &descriptorLayout,
            0U,
            nullptr,
        };
        VkPipelineLayout pipelineLayout = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreatePipelineLayout(device, &pipelineLayoutInfo, nullptr, &pipelineLayout),
            "vkCreatePipelineLayout(counter RNG KAT)");
        pipelineLayout_.reset(
            pipelineLayout,
            [device](const VkPipelineLayout value) {
                vkDestroyPipelineLayout(device, value, nullptr);
            });

        const VkShaderModuleCreateInfo moduleInfo{
            VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,
            nullptr,
            0U,
            shaderCode.size() * sizeof(std::uint32_t),
            shaderCode.data(),
        };
        VkShaderModule shaderModule = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateShaderModule(device, &moduleInfo, nullptr, &shaderModule),
            "vkCreateShaderModule(counter RNG KAT)");
        shaderModule_.reset(
            shaderModule,
            [device](const VkShaderModule value) {
                vkDestroyShaderModule(device, value, nullptr);
            });

        const VkPipelineShaderStageCreateInfo stageInfo{
            VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
            nullptr,
            0U,
            VK_SHADER_STAGE_COMPUTE_BIT,
            shaderModule,
            "main",
            nullptr,
        };
        const VkComputePipelineCreateInfo pipelineInfo{
            VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO,
            nullptr,
            0U,
            stageInfo,
            pipelineLayout,
            VK_NULL_HANDLE,
            -1,
        };
        VkPipeline pipeline = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateComputePipelines(
                device, VK_NULL_HANDLE, 1U, &pipelineInfo, nullptr, &pipeline),
            "vkCreateComputePipelines(counter RNG KAT)");
        pipeline_.reset(
            pipeline,
            [device](const VkPipeline value) {
                vkDestroyPipeline(device, value, nullptr);
            });
    }

    void createDescriptorSet() {
        const VkDevice device = runtime_.device();
        const VkDescriptorPoolSize poolSize{
            VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            1U,
        };
        const VkDescriptorPoolCreateInfo poolInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
            nullptr,
            0U,
            1U,
            1U,
            &poolSize,
        };
        VkDescriptorPool descriptorPool = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateDescriptorPool(device, &poolInfo, nullptr, &descriptorPool),
            "vkCreateDescriptorPool(counter RNG KAT)");
        descriptorPool_.reset(
            descriptorPool,
            [device](const VkDescriptorPool value) {
                vkDestroyDescriptorPool(device, value, nullptr);
            });

        const VkDescriptorSetLayout descriptorLayout = descriptorLayout_.get();
        const VkDescriptorSetAllocateInfo allocateInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
            nullptr,
            descriptorPool,
            1U,
            &descriptorLayout,
        };
        detail::requireVulkanSuccess(
            vkAllocateDescriptorSets(device, &allocateInfo, &descriptorSet_),
            "vkAllocateDescriptorSets(counter RNG KAT)");

        const VkDescriptorBufferInfo outputInfo{
            outputBuffer_.get(),
            0U,
            outputByteSize(),
        };
        const VkWriteDescriptorSet write{
            VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
            nullptr,
            descriptorSet_,
            0U,
            0U,
            1U,
            VK_DESCRIPTOR_TYPE_STORAGE_BUFFER,
            nullptr,
            &outputInfo,
            nullptr,
        };
        vkUpdateDescriptorSets(device, 1U, &write, 0U, nullptr);
    }

    [[nodiscard]] ObservationArray dispatch(const std::uint32_t poisonSalt) {
        std::lock_guard<std::mutex> lock{runtime_.operationMutex()};
        ObservationArray result = makePoisonedObservations(poisonSalt);
        outputBuffer_.write(result.data(), outputByteSize());

        const VkCommandBuffer commands = runtime_.beginCommands();
        const VkMemoryBarrier hostWriteBarrier{
            VK_STRUCTURE_TYPE_MEMORY_BARRIER,
            nullptr,
            VK_ACCESS_HOST_WRITE_BIT,
            VK_ACCESS_SHADER_WRITE_BIT,
        };
        vkCmdPipelineBarrier(
            commands,
            VK_PIPELINE_STAGE_HOST_BIT,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            0U,
            1U,
            &hostWriteBarrier,
            0U,
            nullptr,
            0U,
            nullptr);
        vkCmdBindPipeline(commands, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_.get());
        vkCmdBindDescriptorSets(
            commands,
            VK_PIPELINE_BIND_POINT_COMPUTE,
            pipelineLayout_.get(),
            0U,
            1U,
            &descriptorSet_,
            0U,
            nullptr);
        vkCmdDispatch(commands, kDispatchGroupCount, 1U, 1U);
        const VkMemoryBarrier hostReadBarrier{
            VK_STRUCTURE_TYPE_MEMORY_BARRIER,
            nullptr,
            VK_ACCESS_SHADER_WRITE_BIT,
            VK_ACCESS_HOST_READ_BIT,
        };
        vkCmdPipelineBarrier(
            commands,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            VK_PIPELINE_STAGE_HOST_BIT,
            0U,
            1U,
            &hostReadBarrier,
            0U,
            nullptr,
            0U,
            nullptr);
        runtime_.submitAndWait("Vulkan counter RNG KAT dispatch");
        outputBuffer_.read(result.data(), outputByteSize());
        return result;
    }

    detail::VulkanRuntime runtime_;
    detail::VulkanBuffer outputBuffer_;
    detail::UniqueVulkanHandle<VkDescriptorSetLayout> descriptorLayout_;
    detail::UniqueVulkanHandle<VkPipelineLayout> pipelineLayout_;
    detail::UniqueVulkanHandle<VkShaderModule> shaderModule_;
    detail::UniqueVulkanHandle<VkPipeline> pipeline_;
    detail::UniqueVulkanHandle<VkDescriptorPool> descriptorPool_;
    VkDescriptorSet descriptorSet_{VK_NULL_HANDLE};
};

}  // namespace

VulkanCounterRngKatResult runVulkanCounterRngKat(
    const std::filesystem::path& spirvPath,
    const VulkanRayQueryOptions options) {
    return VulkanCounterRngKatRunner{spirvPath, options}.run();
}

}  // namespace raym0nade::gpu
