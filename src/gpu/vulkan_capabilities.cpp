#include "raym0nade/gpu/vulkan_capabilities.hpp"

#include <vulkan/vulkan.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace raym0nade::gpu {
namespace {

constexpr std::uint32_t kAmdVendorId = 0x1002U;

void requireSuccess(VkResult result, const char* operation) {
    if (result != VK_SUCCESS) {
        throw std::runtime_error(
            std::string{operation} + " failed with Vulkan result " +
            std::to_string(static_cast<int>(result)) + '.');
    }
}

bool versionAtLeast(std::uint32_t version, std::uint32_t major, std::uint32_t minor) noexcept {
    return VK_API_VERSION_MAJOR(version) > major ||
           (VK_API_VERSION_MAJOR(version) == major && VK_API_VERSION_MINOR(version) >= minor);
}

VulkanVersion unpackVersion(std::uint32_t version) noexcept {
    return VulkanVersion{
        VK_API_VERSION_MAJOR(version),
        VK_API_VERSION_MINOR(version),
        VK_API_VERSION_PATCH(version),
    };
}

bool hasExtension(
    const std::vector<VkExtensionProperties>& extensions, const char* requiredName) noexcept {
    return std::any_of(
        extensions.begin(), extensions.end(), [requiredName](const VkExtensionProperties& value) {
            return std::strcmp(value.extensionName, requiredName) == 0;
        });
}

std::vector<VkExtensionProperties> enumerateDeviceExtensions(VkPhysicalDevice device) {
    for (;;) {
        std::uint32_t count = 0;
        requireSuccess(
            vkEnumerateDeviceExtensionProperties(device, nullptr, &count, nullptr),
            "vkEnumerateDeviceExtensionProperties(count)");
        std::vector<VkExtensionProperties> extensions(count);
        if (count == 0) {
            return extensions;
        }
        const VkResult result =
            vkEnumerateDeviceExtensionProperties(device, nullptr, &count, extensions.data());
        if (result == VK_SUCCESS) {
            extensions.resize(count);
            return extensions;
        }
        if (result != VK_INCOMPLETE) {
            requireSuccess(result, "vkEnumerateDeviceExtensionProperties(data)");
        }
    }
}

std::vector<VkExtensionProperties> enumerateInstanceExtensions() {
    for (;;) {
        std::uint32_t count = 0;
        requireSuccess(
            vkEnumerateInstanceExtensionProperties(nullptr, &count, nullptr),
            "vkEnumerateInstanceExtensionProperties(count)");
        std::vector<VkExtensionProperties> extensions(count);
        if (count == 0) {
            return extensions;
        }
        const VkResult result =
            vkEnumerateInstanceExtensionProperties(nullptr, &count, extensions.data());
        if (result == VK_SUCCESS) {
            extensions.resize(count);
            return extensions;
        }
        if (result != VK_INCOMPLETE) {
            requireSuccess(result, "vkEnumerateInstanceExtensionProperties(data)");
        }
    }
}

std::uint64_t deviceLocalMemory(const VkPhysicalDeviceMemoryProperties& properties) noexcept {
    std::uint64_t total = 0;
    for (std::uint32_t index = 0; index < properties.memoryHeapCount; ++index) {
        const VkMemoryHeap& heap = properties.memoryHeaps[index];
        if ((heap.flags & VK_MEMORY_HEAP_DEVICE_LOCAL_BIT) == 0) {
            continue;
        }
        if (heap.size > std::numeric_limits<std::uint64_t>::max() - total) {
            return std::numeric_limits<std::uint64_t>::max();
        }
        total += heap.size;
    }
    return total;
}

std::pair<bool, std::uint32_t> findComputeQueue(VkPhysicalDevice device) {
    std::uint32_t count = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(device, &count, nullptr);
    std::vector<VkQueueFamilyProperties> families(count);
    if (count != 0) {
        vkGetPhysicalDeviceQueueFamilyProperties(device, &count, families.data());
        families.resize(count);
    }
    for (std::uint32_t index = 0; index < families.size(); ++index) {
        if (families[index].queueCount != 0 &&
            (families[index].queueFlags & VK_QUEUE_COMPUTE_BIT) != 0) {
            return {true, index};
        }
    }
    return {false, 0};
}

void appendMissingRequirements(
    VulkanDeviceCapabilities& result, bool loaderSupportsVulkan12) {
    if (!loaderSupportsVulkan12) {
        result.missingRequirements.emplace_back("Vulkan 1.2 loader");
    }
    if (!versionAtLeast(
            VK_MAKE_API_VERSION(
                0, result.apiVersion.major, result.apiVersion.minor, result.apiVersion.patch),
            1,
            2)) {
        result.missingRequirements.emplace_back("Vulkan 1.2 device");
    }
    if (!result.hasComputeQueue) {
        result.missingRequirements.emplace_back("compute queue");
    }
    if (!result.bufferDeviceAddress) {
        result.missingRequirements.emplace_back("buffer device address");
    }
    if (!result.accelerationStructure) {
        result.missingRequirements.emplace_back("acceleration structure");
    }
    if (!result.rayQuery) {
        result.missingRequirements.emplace_back("ray query");
    }
}

VulkanDeviceCapabilities inspectDevice(VkPhysicalDevice device, bool loaderSupportsVulkan12) {
    const std::vector<VkExtensionProperties> extensions = enumerateDeviceExtensions(device);
    const bool hasAccelerationStructureExtension =
        hasExtension(extensions, VK_KHR_ACCELERATION_STRUCTURE_EXTENSION_NAME);
    const bool hasDeferredHostOperationsExtension =
        hasExtension(extensions, VK_KHR_DEFERRED_HOST_OPERATIONS_EXTENSION_NAME);
    const bool hasRayQueryExtension = hasExtension(extensions, VK_KHR_RAY_QUERY_EXTENSION_NAME);
    const bool hasRayTracingPipelineExtension =
        hasExtension(extensions, VK_KHR_RAY_TRACING_PIPELINE_EXTENSION_NAME);

    VkPhysicalDeviceProperties basicProperties{};
    vkGetPhysicalDeviceProperties(device, &basicProperties);

    VkPhysicalDeviceMemoryProperties memoryProperties{};
    vkGetPhysicalDeviceMemoryProperties(device, &memoryProperties);
    const auto [hasComputeQueue, computeQueueFamily] = findComputeQueue(device);

    VulkanDeviceCapabilities result;
    result.deviceName = basicProperties.deviceName;
    result.apiVersion = unpackVersion(basicProperties.apiVersion);
    result.deviceLocalMemoryBytes = deviceLocalMemory(memoryProperties);
    result.vendorId = basicProperties.vendorID;
    result.deviceId = basicProperties.deviceID;
    result.computeQueueFamily = computeQueueFamily;
    result.hasComputeQueue = hasComputeQueue;
    result.integrated = basicProperties.deviceType == VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU;

    const bool deviceSupportsVulkan12 = versionAtLeast(basicProperties.apiVersion, 1, 2);
    if (!loaderSupportsVulkan12 || !deviceSupportsVulkan12) {
        result.driverName = "unavailable";
        result.driverInfo = "Vulkan 1.2 properties unavailable";
        appendMissingRequirements(result, loaderSupportsVulkan12);
        return result;
    }

    VkPhysicalDeviceRayTracingPipelineFeaturesKHR rayTracingPipelineFeatures{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_RAY_TRACING_PIPELINE_FEATURES_KHR};
    VkPhysicalDeviceRayQueryFeaturesKHR rayQueryFeatures{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_RAY_QUERY_FEATURES_KHR};
    VkPhysicalDeviceAccelerationStructureFeaturesKHR accelerationStructureFeatures{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ACCELERATION_STRUCTURE_FEATURES_KHR};
    VkPhysicalDeviceBufferDeviceAddressFeatures bufferDeviceAddressFeatures{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_BUFFER_DEVICE_ADDRESS_FEATURES};
    VkPhysicalDeviceFeatures2 features{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
    features.pNext = &bufferDeviceAddressFeatures;
    void** featureTail = &bufferDeviceAddressFeatures.pNext;
    if (hasAccelerationStructureExtension) {
        *featureTail = &accelerationStructureFeatures;
        featureTail = &accelerationStructureFeatures.pNext;
    }
    if (hasRayQueryExtension) {
        *featureTail = &rayQueryFeatures;
        featureTail = &rayQueryFeatures.pNext;
    }
    if (hasRayTracingPipelineExtension) {
        *featureTail = &rayTracingPipelineFeatures;
    }
    vkGetPhysicalDeviceFeatures2(device, &features);

    VkPhysicalDeviceAccelerationStructurePropertiesKHR accelerationStructureProperties{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ACCELERATION_STRUCTURE_PROPERTIES_KHR};
    VkPhysicalDeviceSubgroupProperties subgroupProperties{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_PROPERTIES};
    VkPhysicalDeviceDriverProperties driverProperties{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_DRIVER_PROPERTIES};
    VkPhysicalDeviceProperties2 properties{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
    properties.pNext = &driverProperties;
    driverProperties.pNext = &subgroupProperties;
    if (hasAccelerationStructureExtension) {
        subgroupProperties.pNext = &accelerationStructureProperties;
    }
    vkGetPhysicalDeviceProperties2(device, &properties);

    result.driverName = driverProperties.driverName;
    result.driverInfo = driverProperties.driverInfo;
    result.maximumAccelerationStructurePrimitiveCount =
        accelerationStructureProperties.maxPrimitiveCount;
    result.subgroupSize = subgroupProperties.subgroupSize;
    result.bufferDeviceAddress = bufferDeviceAddressFeatures.bufferDeviceAddress == VK_TRUE;
    result.accelerationStructure = hasAccelerationStructureExtension &&
                                   hasDeferredHostOperationsExtension &&
                                   accelerationStructureFeatures.accelerationStructure == VK_TRUE;
    result.rayQuery = hasRayQueryExtension && rayQueryFeatures.rayQuery == VK_TRUE;
    result.rayTracingPipeline = hasRayTracingPipelineExtension &&
                                rayTracingPipelineFeatures.rayTracingPipeline == VK_TRUE;

    appendMissingRequirements(result, loaderSupportsVulkan12);
    return result;
}

class Instance {
public:
    explicit Instance(std::uint32_t apiVersion) {
        const std::vector<VkExtensionProperties> extensions = enumerateInstanceExtensions();
        const bool enumeratePortability =
            hasExtension(extensions, VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME);
        const char* enabledExtensions[] = {VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME};
        const VkApplicationInfo applicationInfo{
            VK_STRUCTURE_TYPE_APPLICATION_INFO,
            nullptr,
            "Raym0nade GPU Probe",
            VK_MAKE_API_VERSION(0, 1, 0, 0),
            "Raym0nade",
            VK_MAKE_API_VERSION(0, 1, 0, 0),
            apiVersion,
        };
        const VkInstanceCreateInfo createInfo{
            VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO,
            nullptr,
            enumeratePortability ? VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR : 0U,
            &applicationInfo,
            0,
            nullptr,
            enumeratePortability ? 1U : 0U,
            enumeratePortability ? enabledExtensions : nullptr,
        };
        requireSuccess(vkCreateInstance(&createInfo, nullptr, &handle_), "vkCreateInstance");
    }

    Instance(const Instance&) = delete;
    Instance& operator=(const Instance&) = delete;

    ~Instance() {
        if (handle_ != VK_NULL_HANDLE) {
            vkDestroyInstance(handle_, nullptr);
        }
    }

    [[nodiscard]] VkInstance get() const noexcept {
        return handle_;
    }

private:
    VkInstance handle_{VK_NULL_HANDLE};
};

}  // namespace

bool VulkanDeviceCapabilities::isAmd() const noexcept {
    return vendorId == kAmdVendorId;
}

bool VulkanDeviceCapabilities::supportsRayQueryBackend() const noexcept {
    const bool supportsVulkan12 =
        apiVersion.major > 1 || (apiVersion.major == 1 && apiVersion.minor >= 2);
    return isAmd() && supportsVulkan12 && hasComputeQueue && bufferDeviceAddress &&
           accelerationStructure && rayQuery;
}

VulkanVersion vulkanLoaderVersion() {
    std::uint32_t version = VK_API_VERSION_1_0;
    const auto enumerateVersion = reinterpret_cast<PFN_vkEnumerateInstanceVersion>(
        vkGetInstanceProcAddr(VK_NULL_HANDLE, "vkEnumerateInstanceVersion"));
    if (enumerateVersion != nullptr) {
        requireSuccess(enumerateVersion(&version), "vkEnumerateInstanceVersion");
    }
    return unpackVersion(version);
}

std::vector<VulkanDeviceCapabilities> enumerateVulkanDevices() {
    const VulkanVersion loader = vulkanLoaderVersion();
    const bool loaderSupportsVulkan12 =
        loader.major > 1 || (loader.major == 1 && loader.minor >= 2);
    const std::uint32_t requestedVersion = loaderSupportsVulkan12
                                               ? VK_API_VERSION_1_2
                                               : VK_MAKE_API_VERSION(0, loader.major, loader.minor, 0);
    Instance instance{requestedVersion};

    std::vector<VkPhysicalDevice> physicalDevices;
    for (;;) {
        std::uint32_t count = 0;
        requireSuccess(
            vkEnumeratePhysicalDevices(instance.get(), &count, nullptr),
            "vkEnumeratePhysicalDevices(count)");
        physicalDevices.resize(count);
        if (count == 0) {
            break;
        }
        const VkResult enumerationResult =
            vkEnumeratePhysicalDevices(instance.get(), &count, physicalDevices.data());
        if (enumerationResult == VK_SUCCESS) {
            physicalDevices.resize(count);
            break;
        }
        if (enumerationResult != VK_INCOMPLETE) {
            requireSuccess(enumerationResult, "vkEnumeratePhysicalDevices(data)");
        }
    }

    std::vector<VulkanDeviceCapabilities> result;
    result.reserve(physicalDevices.size());
    for (VkPhysicalDevice device : physicalDevices) {
        result.push_back(inspectDevice(device, loaderSupportsVulkan12));
    }
    return result;
}

}  // namespace raym0nade::gpu
