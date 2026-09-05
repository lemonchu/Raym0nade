#include "raym0nade/gpu/vulkan_ray_query.hpp"

#include <vulkan/vulkan.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include "raym0nade/scene_data.hpp"
#include "vulkan_runtime.hpp"

namespace raym0nade::gpu {
namespace {

constexpr std::uint32_t kWorkgroupSize = 64U;

bool allFinite(const std::array<float, 4>& values) noexcept {
    return std::all_of(values.begin(), values.end(), [](float value) {
        return std::isfinite(value);
    });
}

void validateRays(const std::vector<VulkanRayQueryRay>& rays) {
    for (const VulkanRayQueryRay& ray : rays) {
        if (!allFinite(ray.originAndTMin) || !allFinite(ray.directionAndTMax)) {
            throw std::invalid_argument("Vulkan Ray Query inputs must contain finite values.");
        }
        const float tMin = ray.originAndTMin[3];
        const float tMax = ray.directionAndTMax[3];
        if (tMin < 0.0F || tMax < tMin) {
            throw std::invalid_argument("Vulkan Ray Query inputs have an invalid ray interval.");
        }
        if (ray.directionAndTMax[0] == 0.0F && ray.directionAndTMax[1] == 0.0F &&
            ray.directionAndTMax[2] == 0.0F) {
            throw std::invalid_argument("Vulkan Ray Query directions must be non-zero.");
        }
    }
}

const PackedSceneData& rejectCutoutIntersectorScene(const PackedSceneData& scene) {
    for (std::uint32_t materialId : scene.triangleMaterialIds) {
        if (materialId < scene.materials.size() &&
            (scene.materials[materialId].flagsAndReserved[0] &
             kPackedMaterialCutout) != 0U) {
            throw std::invalid_argument(
                "VulkanRayQueryIntersector does not support alpha-cutout materials.");
        }
    }
    return scene;
}

}  // namespace

class VulkanRayQueryIntersector::Implementation {
public:
    Implementation(
        const PackedSceneData& scene,
        const std::filesystem::path& spirvPath,
        VulkanRayQueryOptions options)
        : runtime_(rejectCutoutIntersectorScene(scene), options) {
        createPipeline(detail::readSpirvFile(spirvPath));
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

    [[nodiscard]] VulkanRayQueryBatch intersect(
        const std::vector<VulkanRayQueryRay>& rays) {
        std::lock_guard<std::mutex> lock{runtime_.operationMutex()};
        VulkanRayQueryBatch result;
        if (rays.empty()) {
            return result;
        }
        validateRays(rays);
        validateBatchSize(rays.size());
        ensureBatchCapacity(rays.size());

        const detail::VulkanClock::time_point begin = detail::VulkanClock::now();
        const VkDeviceSize rayBytes = detail::checkedVulkanByteSize(
            rays.size(), sizeof(VulkanRayQueryRay), "Vulkan ray batch");
        rayBuffer_->write(rays.data(), rayBytes);

        const VkCommandBuffer commands = runtime_.beginCommands();
        const VkMemoryBarrier uploadBarrier{
            VK_STRUCTURE_TYPE_MEMORY_BARRIER,
            nullptr,
            VK_ACCESS_HOST_WRITE_BIT,
            VK_ACCESS_SHADER_READ_BIT,
        };
        vkCmdPipelineBarrier(
            commands,
            VK_PIPELINE_STAGE_HOST_BIT,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            0U,
            1U,
            &uploadBarrier,
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
        const std::uint32_t rayCount = static_cast<std::uint32_t>(rays.size());
        vkCmdPushConstants(
            commands,
            pipelineLayout_.get(),
            VK_SHADER_STAGE_COMPUTE_BIT,
            0U,
            sizeof(rayCount),
            &rayCount);
        const std::uint32_t groupCount = static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(rayCount) + kWorkgroupSize - 1U) /
            kWorkgroupSize);
        vkCmdDispatch(commands, groupCount, 1U, 1U);
        const VkMemoryBarrier readbackBarrier{
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
            &readbackBarrier,
            0U,
            nullptr,
            0U,
            nullptr);
        runtime_.submitAndWait("Vulkan Ray Query batch");

        result.hits.resize(rays.size());
        const VkDeviceSize resultBytes = detail::checkedVulkanByteSize(
            result.hits.size(),
            sizeof(VulkanRayQueryHit),
            "Vulkan Ray Query result batch");
        hitBuffer_->read(result.hits.data(), resultBytes);
        result.dispatchAndReadbackMilliseconds = detail::elapsedMilliseconds(begin);
        validateHits(result.hits);
        return result;
    }

private:
    void createPipeline(const std::vector<std::uint32_t>& shaderCode) {
        const VkDevice device = runtime_.device();
        std::array<VkDescriptorSetLayoutBinding, 3> bindings{};
        bindings[0].binding = 0U;
        bindings[0].descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
        bindings[0].descriptorCount = 1U;
        bindings[0].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        bindings[1].binding = 1U;
        bindings[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        bindings[1].descriptorCount = 1U;
        bindings[1].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        bindings[2].binding = 2U;
        bindings[2].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        bindings[2].descriptorCount = 1U;
        bindings[2].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        const VkDescriptorSetLayoutCreateInfo layoutInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,
            nullptr,
            0U,
            static_cast<std::uint32_t>(bindings.size()),
            bindings.data(),
        };
        VkDescriptorSetLayout descriptorLayout = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateDescriptorSetLayout(device, &layoutInfo, nullptr, &descriptorLayout),
            "vkCreateDescriptorSetLayout");
        descriptorLayout_.reset(
            descriptorLayout,
            [device](VkDescriptorSetLayout value) {
                vkDestroyDescriptorSetLayout(device, value, nullptr);
            });

        const VkPushConstantRange pushConstant{
            VK_SHADER_STAGE_COMPUTE_BIT,
            0U,
            sizeof(std::uint32_t),
        };
        const VkPipelineLayoutCreateInfo pipelineLayoutInfo{
            VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,
            nullptr,
            0U,
            1U,
            &descriptorLayout,
            1U,
            &pushConstant,
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
            "vkCreateShaderModule");
        shaderModule_.reset(
            shaderModule,
            [device](VkShaderModule value) {
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
            "vkCreateComputePipelines");
        pipeline_.reset(
            pipeline,
            [device](VkPipeline value) { vkDestroyPipeline(device, value, nullptr); });

        const std::array<VkDescriptorPoolSize, 2> poolSizes{{
            {VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR, 1U},
            {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 2U},
        }};
        const VkDescriptorPoolCreateInfo poolInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
            nullptr,
            0U,
            1U,
            static_cast<std::uint32_t>(poolSizes.size()),
            poolSizes.data(),
        };
        VkDescriptorPool descriptorPool = VK_NULL_HANDLE;
        detail::requireVulkanSuccess(
            vkCreateDescriptorPool(device, &poolInfo, nullptr, &descriptorPool),
            "vkCreateDescriptorPool");
        descriptorPool_.reset(
            descriptorPool,
            [device](VkDescriptorPool value) {
                vkDestroyDescriptorPool(device, value, nullptr);
            });
        const VkDescriptorSetAllocateInfo allocateInfo{
            VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
            nullptr,
            descriptorPool,
            1U,
            &descriptorLayout,
        };
        detail::requireVulkanSuccess(
            vkAllocateDescriptorSets(device, &allocateInfo, &descriptorSet_),
            "vkAllocateDescriptorSets");
    }

    void validateBatchSize(std::size_t rayCount) const {
        if (rayCount > std::numeric_limits<std::uint32_t>::max()) {
            throw std::overflow_error(
                "A Vulkan Ray Query batch exceeds the 32-bit shader ABI.");
        }
        const VkDeviceSize rayBytes = detail::checkedVulkanByteSize(
            rayCount, sizeof(VulkanRayQueryRay), "Vulkan ray batch");
        const VkDeviceSize hitBytes = detail::checkedVulkanByteSize(
            rayCount, sizeof(VulkanRayQueryHit), "Vulkan hit batch");
        const VkPhysicalDeviceLimits& limits = runtime_.physicalProperties().limits;
        if (rayBytes > limits.maxStorageBufferRange ||
            hitBytes > limits.maxStorageBufferRange) {
            throw std::invalid_argument(
                "The Vulkan Ray Query batch exceeds maxStorageBufferRange.");
        }
        const std::uint64_t groups =
            (static_cast<std::uint64_t>(rayCount) + kWorkgroupSize - 1U) /
            kWorkgroupSize;
        if (groups > limits.maxComputeWorkGroupCount[0]) {
            throw std::invalid_argument(
                "The Vulkan Ray Query batch exceeds the X dispatch limit.");
        }
    }

    void ensureBatchCapacity(std::size_t rayCount) {
        if (rayCount <= batchCapacity_) {
            return;
        }
        const VkDeviceSize rayBytes = detail::checkedVulkanByteSize(
            rayCount, sizeof(VulkanRayQueryRay), "Vulkan ray batch");
        const VkDeviceSize hitBytes = detail::checkedVulkanByteSize(
            rayCount, sizeof(VulkanRayQueryHit), "Vulkan hit batch");
        constexpr VkMemoryPropertyFlags kHostMemory =
            VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
            VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
        const VkDevice device = runtime_.device();
        auto newRayBuffer = std::make_unique<detail::VulkanBuffer>(
            device,
            runtime_.memoryProperties(),
            rayBytes,
            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
            kHostMemory,
            false);
        auto newHitBuffer = std::make_unique<detail::VulkanBuffer>(
            device,
            runtime_.memoryProperties(),
            hitBytes,
            VK_BUFFER_USAGE_STORAGE_BUFFER_BIT,
            kHostMemory,
            false);

        const VkAccelerationStructureKHR topLevelHandle =
            runtime_.topLevelAccelerationStructure();
        const VkWriteDescriptorSetAccelerationStructureKHR accelerationWrite{
            VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET_ACCELERATION_STRUCTURE_KHR,
            nullptr,
            1U,
            &topLevelHandle,
        };
        const VkDescriptorBufferInfo rayInfo{newRayBuffer->get(), 0U, rayBytes};
        const VkDescriptorBufferInfo hitInfo{newHitBuffer->get(), 0U, hitBytes};
        std::array<VkWriteDescriptorSet, 3> writes{};
        writes[0].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[0].pNext = &accelerationWrite;
        writes[0].dstSet = descriptorSet_;
        writes[0].dstBinding = 0U;
        writes[0].descriptorCount = 1U;
        writes[0].descriptorType = VK_DESCRIPTOR_TYPE_ACCELERATION_STRUCTURE_KHR;
        writes[1].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[1].dstSet = descriptorSet_;
        writes[1].dstBinding = 1U;
        writes[1].descriptorCount = 1U;
        writes[1].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[1].pBufferInfo = &rayInfo;
        writes[2].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[2].dstSet = descriptorSet_;
        writes[2].dstBinding = 2U;
        writes[2].descriptorCount = 1U;
        writes[2].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[2].pBufferInfo = &hitInfo;
        vkUpdateDescriptorSets(
            device,
            static_cast<std::uint32_t>(writes.size()),
            writes.data(),
            0U,
            nullptr);
        rayBuffer_ = std::move(newRayBuffer);
        hitBuffer_ = std::move(newHitBuffer);
        batchCapacity_ = rayCount;
    }

    void validateHits(const std::vector<VulkanRayQueryHit>& hits) const {
        for (const VulkanRayQueryHit& hit : hits) {
            if (hit.hit > 1U) {
                throw std::runtime_error("The Vulkan shader returned an invalid hit flag.");
            }
            if (hit.hasHit()) {
                if (hit.primitiveId >= runtime_.triangleCount() ||
                    !std::isfinite(hit.distance) ||
                    !std::isfinite(hit.barycentricU) ||
                    !std::isfinite(hit.barycentricV)) {
                    throw std::runtime_error(
                        "The Vulkan shader returned an invalid hit record.");
                }
            } else if (hit.primitiveId != kInvalidRayQueryPrimitiveId) {
                throw std::runtime_error(
                    "The Vulkan shader returned an invalid miss sentinel.");
            }
        }
    }

    detail::VulkanRuntime runtime_;
    std::size_t batchCapacity_{0U};
    std::unique_ptr<detail::VulkanBuffer> rayBuffer_;
    std::unique_ptr<detail::VulkanBuffer> hitBuffer_;

    detail::UniqueVulkanHandle<VkDescriptorSetLayout> descriptorLayout_;
    detail::UniqueVulkanHandle<VkPipelineLayout> pipelineLayout_;
    detail::UniqueVulkanHandle<VkShaderModule> shaderModule_;
    detail::UniqueVulkanHandle<VkPipeline> pipeline_;
    detail::UniqueVulkanHandle<VkDescriptorPool> descriptorPool_;
    VkDescriptorSet descriptorSet_{VK_NULL_HANDLE};
};

VulkanRayQueryIntersector::VulkanRayQueryIntersector(
    const PackedSceneData& scene,
    const std::filesystem::path& spirvPath,
    VulkanRayQueryOptions options)
    : implementation_(std::make_unique<Implementation>(scene, spirvPath, options)) {}

VulkanRayQueryIntersector::~VulkanRayQueryIntersector() = default;

const std::string& VulkanRayQueryIntersector::deviceName() const noexcept {
    return implementation_->deviceName();
}

const VulkanRayQuerySetupTimings&
VulkanRayQueryIntersector::setupTimings() const noexcept {
    return implementation_->setupTimings();
}

VulkanValidationReport VulkanRayQueryIntersector::validationReport() const {
    return implementation_->validationReport();
}

VulkanRayQueryBatch VulkanRayQueryIntersector::intersect(
    const std::vector<VulkanRayQueryRay>& rays) {
    return implementation_->intersect(rays);
}

}  // namespace raym0nade::gpu
