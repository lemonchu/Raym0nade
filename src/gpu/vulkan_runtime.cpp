#include "vulkan_runtime.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "raym0nade/scene_data.hpp"

namespace raym0nade::gpu::detail {
namespace {

constexpr std::uint32_t kAmdVendorId = 0x1002U;
constexpr std::uint32_t kInstanceMask = 0xffU;
constexpr VkDeviceSize kPackedUploadChunkBytes = 16U * 1024U * 1024U;

bool versionAtLeast(
    std::uint32_t version, std::uint32_t major, std::uint32_t minor) noexcept {
    return VK_API_VERSION_MAJOR(version) > major ||
           (VK_API_VERSION_MAJOR(version) == major &&
            VK_API_VERSION_MINOR(version) >= minor);
}

template <typename Property, typename Enumerate>
std::vector<Property> enumerateProperties(Enumerate&& enumerate, const char* operation) {
    for (;;) {
        std::uint32_t count = 0U;
        const VkResult countResult = enumerate(&count, nullptr);
        if (countResult != VK_SUCCESS && countResult != VK_INCOMPLETE) {
            requireVulkanSuccess(countResult, operation);
        }
        std::vector<Property> properties(count);
        if (count == 0U) {
            if (countResult == VK_INCOMPLETE) {
                continue;
            }
            return properties;
        }
        const VkResult result = enumerate(&count, properties.data());
        if (result == VK_SUCCESS) {
            properties.resize(count);
            return properties;
        }
        if (result != VK_INCOMPLETE) {
            requireVulkanSuccess(result, operation);
        }
    }
}

std::vector<VkLayerProperties> enumerateInstanceLayers() {
    return enumerateProperties<VkLayerProperties>(
        [](std::uint32_t* count, VkLayerProperties* values) {
            return vkEnumerateInstanceLayerProperties(count, values);
        },
        "vkEnumerateInstanceLayerProperties");
}

std::vector<VkExtensionProperties> enumerateInstanceExtensions(
    const char* layerName = nullptr) {
    return enumerateProperties<VkExtensionProperties>(
        [layerName](std::uint32_t* count, VkExtensionProperties* values) {
            return vkEnumerateInstanceExtensionProperties(layerName, count, values);
        },
        "vkEnumerateInstanceExtensionProperties");
}

std::vector<VkExtensionProperties> enumerateDeviceExtensions(
    VkPhysicalDevice physicalDevice) {
    return enumerateProperties<VkExtensionProperties>(
        [physicalDevice](std::uint32_t* count, VkExtensionProperties* values) {
            return vkEnumerateDeviceExtensionProperties(
                physicalDevice, nullptr, count, values);
        },
        "vkEnumerateDeviceExtensionProperties");
}

bool hasLayer(const std::vector<VkLayerProperties>& layers, const char* name) noexcept {
    return std::any_of(layers.begin(), layers.end(), [name](const VkLayerProperties& layer) {
        return std::strcmp(layer.layerName, name) == 0;
    });
}

bool hasExtension(
    const std::vector<VkExtensionProperties>& extensions, const char* name) noexcept {
    return std::any_of(
        extensions.begin(),
        extensions.end(),
        [name](const VkExtensionProperties& extension) {
            return std::strcmp(extension.extensionName, name) == 0;
        });
}

struct ValidationState {
    explicit ValidationState(bool wasRequested) : requested(wasRequested) {}

    std::atomic<std::uint32_t> errorCount{0U};
    std::atomic<std::uint32_t> warningCount{0U};
    mutable std::mutex messageMutex;
    std::vector<std::string> messages;
    bool requested{false};
    bool enabled{false};
    bool synchronizationEnabled{false};

    void record(
        VkDebugUtilsMessageSeverityFlagBitsEXT severity, const char* message) noexcept {
        if ((severity & VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT) != 0U) {
            errorCount.fetch_add(1U, std::memory_order_relaxed);
        } else if ((severity & VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT) != 0U) {
            warningCount.fetch_add(1U, std::memory_order_relaxed);
        }
        if (message == nullptr) {
            return;
        }
        try {
            std::lock_guard<std::mutex> lock{messageMutex};
            constexpr std::size_t kMaximumStoredMessages = 64U;
            if (messages.size() < kMaximumStoredMessages) {
                messages.emplace_back(message);
            }
        } catch (...) {
            // A Vulkan callback must never propagate an exception across the C ABI boundary.
        }
    }

    [[nodiscard]] VulkanValidationReport snapshot() const {
        VulkanValidationReport result;
        result.requested = requested;
        result.enabled = enabled;
        result.synchronizationValidationEnabled = synchronizationEnabled;
        result.errorCount = errorCount.load(std::memory_order_relaxed);
        result.warningCount = warningCount.load(std::memory_order_relaxed);
        std::lock_guard<std::mutex> lock{messageMutex};
        result.messages = messages;
        return result;
    }
};

VKAPI_ATTR VkBool32 VKAPI_CALL validationCallback(
    VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity,
    VkDebugUtilsMessageTypeFlagsEXT,
    const VkDebugUtilsMessengerCallbackDataEXT* callbackData,
    void* userData) noexcept {
    auto* state = static_cast<ValidationState*>(userData);
    if (state != nullptr) {
        state->record(
            messageSeverity, callbackData == nullptr ? nullptr : callbackData->pMessage);
    }
    return VK_FALSE;
}

class VulkanInstance {
public:
    VulkanInstance(bool requestValidation, ValidationState& validation) {
        std::uint32_t loaderVersion = VK_API_VERSION_1_0;
        const auto enumerateVersion = reinterpret_cast<PFN_vkEnumerateInstanceVersion>(
            vkGetInstanceProcAddr(VK_NULL_HANDLE, "vkEnumerateInstanceVersion"));
        if (enumerateVersion != nullptr) {
            requireVulkanSuccess(
                enumerateVersion(&loaderVersion), "vkEnumerateInstanceVersion");
        }
        if (!versionAtLeast(loaderVersion, 1U, 2U)) {
            throw std::runtime_error("The Vulkan loader does not support Vulkan 1.2.");
        }

        const std::vector<VkLayerProperties> layers = enumerateInstanceLayers();
        const std::vector<VkExtensionProperties> extensions = enumerateInstanceExtensions();
        const bool hasValidationLayer =
            hasLayer(layers, "VK_LAYER_KHRONOS_validation");
        const bool hasDebugUtils =
            hasExtension(extensions, VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
        const std::vector<VkExtensionProperties> validationLayerExtensions =
            hasValidationLayer
                ? enumerateInstanceExtensions("VK_LAYER_KHRONOS_validation")
                : std::vector<VkExtensionProperties>{};
        const bool hasValidationFeatures = hasExtension(
            validationLayerExtensions, VK_EXT_VALIDATION_FEATURES_EXTENSION_NAME);
        validation.enabled = requestValidation && hasValidationLayer && hasDebugUtils;
        validation.synchronizationEnabled = validation.enabled && hasValidationFeatures;

        std::vector<const char*> enabledLayers;
        std::vector<const char*> enabledExtensions;
        if (validation.enabled) {
            enabledLayers.push_back("VK_LAYER_KHRONOS_validation");
            enabledExtensions.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);
            if (validation.synchronizationEnabled) {
                enabledExtensions.push_back(VK_EXT_VALIDATION_FEATURES_EXTENSION_NAME);
            }
        }

        VkInstanceCreateFlags flags = 0U;
#ifdef VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME
        if (hasExtension(extensions, VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME)) {
            enabledExtensions.push_back(VK_KHR_PORTABILITY_ENUMERATION_EXTENSION_NAME);
            flags |= VK_INSTANCE_CREATE_ENUMERATE_PORTABILITY_BIT_KHR;
        }
#endif

        const VkApplicationInfo applicationInfo{
            VK_STRUCTURE_TYPE_APPLICATION_INFO,
            nullptr,
            "Raym0nade Vulkan Ray Query",
            VK_MAKE_API_VERSION(0, 1, 0, 0),
            "Raym0nade",
            VK_MAKE_API_VERSION(0, 1, 0, 0),
            VK_API_VERSION_1_2,
        };
        VkDebugUtilsMessengerCreateInfoEXT debugCreateInfo{
            VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT};
        debugCreateInfo.messageSeverity =
            VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT |
            VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT;
        debugCreateInfo.messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT |
                                      VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT |
                                      VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT;
        debugCreateInfo.pfnUserCallback = validationCallback;
        debugCreateInfo.pUserData = &validation;

        const VkValidationFeatureEnableEXT synchronizationFeature =
            VK_VALIDATION_FEATURE_ENABLE_SYNCHRONIZATION_VALIDATION_EXT;
        VkValidationFeaturesEXT validationFeatures{
            VK_STRUCTURE_TYPE_VALIDATION_FEATURES_EXT};
        validationFeatures.enabledValidationFeatureCount = 1U;
        validationFeatures.pEnabledValidationFeatures = &synchronizationFeature;
        validationFeatures.pNext = &debugCreateInfo;

        const void* instanceChain = nullptr;
        if (validation.synchronizationEnabled) {
            instanceChain = &validationFeatures;
        } else if (validation.enabled) {
            instanceChain = &debugCreateInfo;
        }
        const VkInstanceCreateInfo createInfo{
            VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO,
            instanceChain,
            flags,
            &applicationInfo,
            static_cast<std::uint32_t>(enabledLayers.size()),
            enabledLayers.empty() ? nullptr : enabledLayers.data(),
            static_cast<std::uint32_t>(enabledExtensions.size()),
            enabledExtensions.empty() ? nullptr : enabledExtensions.data(),
        };

        VkInstance instance = VK_NULL_HANDLE;
        requireVulkanSuccess(
            vkCreateInstance(&createInfo, nullptr, &instance), "vkCreateInstance");
        instance_.reset(instance, [](VkInstance value) { vkDestroyInstance(value, nullptr); });

        if (validation.enabled) {
            const auto createMessenger =
                reinterpret_cast<PFN_vkCreateDebugUtilsMessengerEXT>(
                    vkGetInstanceProcAddr(instance, "vkCreateDebugUtilsMessengerEXT"));
            const auto destroyMessenger =
                reinterpret_cast<PFN_vkDestroyDebugUtilsMessengerEXT>(
                    vkGetInstanceProcAddr(instance, "vkDestroyDebugUtilsMessengerEXT"));
            if (createMessenger == nullptr || destroyMessenger == nullptr) {
                throw std::runtime_error(
                    "The Vulkan debug-utils entry points are unavailable.");
            }
            VkDebugUtilsMessengerEXT messenger = VK_NULL_HANDLE;
            requireVulkanSuccess(
                createMessenger(instance, &debugCreateInfo, nullptr, &messenger),
                "vkCreateDebugUtilsMessengerEXT");
            messenger_.reset(
                messenger,
                [instance, destroyMessenger](VkDebugUtilsMessengerEXT value) {
                    destroyMessenger(instance, value, nullptr);
                });
        }
    }

    VulkanInstance(const VulkanInstance&) = delete;
    VulkanInstance& operator=(const VulkanInstance&) = delete;

    [[nodiscard]] VkInstance get() const noexcept {
        return instance_.get();
    }

private:
    UniqueVulkanHandle<VkInstance> instance_;
    UniqueVulkanHandle<VkDebugUtilsMessengerEXT> messenger_;
};

std::pair<bool, std::uint32_t> findComputeQueue(VkPhysicalDevice physicalDevice) {
    std::uint32_t count = 0U;
    vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, &count, nullptr);
    std::vector<VkQueueFamilyProperties> families(count);
    if (count != 0U) {
        vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, &count, families.data());
        families.resize(count);
    }
    for (std::uint32_t index = 0U; index < count; ++index) {
        const VkQueueFamilyProperties& family = families[index];
        if (family.queueCount != 0U &&
            (family.queueFlags & VK_QUEUE_COMPUTE_BIT) != 0U &&
            (family.queueFlags & VK_QUEUE_GRAPHICS_BIT) == 0U) {
            return {true, index};
        }
    }
    for (std::uint32_t index = 0U; index < count; ++index) {
        const VkQueueFamilyProperties& family = families[index];
        if (family.queueCount != 0U &&
            (family.queueFlags & VK_QUEUE_COMPUTE_BIT) != 0U) {
            return {true, index};
        }
    }
    return {false, 0U};
}

struct PhysicalDeviceSelection {
    VkPhysicalDevice physicalDevice{VK_NULL_HANDLE};
    VkPhysicalDeviceMemoryProperties memoryProperties{};
    VkPhysicalDeviceProperties properties{};
    VkPhysicalDeviceAccelerationStructurePropertiesKHR accelerationProperties{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ACCELERATION_STRUCTURE_PROPERTIES_KHR};
    VkDeviceSize scratchAlignment{1U};
    std::uint32_t queueFamily{0U};
    std::uint32_t timestampValidBits{0U};
    std::string deviceName;
    bool requiresPortabilitySubset{false};
};

std::uint32_t queueTimestampValidBits(
    VkPhysicalDevice physicalDevice, std::uint32_t queueFamily) {
    std::uint32_t count = 0U;
    vkGetPhysicalDeviceQueueFamilyProperties(physicalDevice, &count, nullptr);
    std::vector<VkQueueFamilyProperties> families(count);
    if (count != 0U) {
        vkGetPhysicalDeviceQueueFamilyProperties(
            physicalDevice, &count, families.data());
        families.resize(count);
    }
    if (!(queueFamily < families.size())) {
        throw std::runtime_error("The selected Vulkan queue family is out of range.");
    }
    return families[queueFamily].timestampValidBits;
}

PhysicalDeviceSelection selectAmdDevice(VkInstance instance) {
    const std::vector<VkPhysicalDevice> devices = enumerateProperties<VkPhysicalDevice>(
        [instance](std::uint32_t* count, VkPhysicalDevice* values) {
            return vkEnumeratePhysicalDevices(instance, count, values);
        },
        "vkEnumeratePhysicalDevices");

    PhysicalDeviceSelection selected;
    int selectedScore = -1;
    std::vector<std::string> rejectedAmdDevices;
    for (VkPhysicalDevice physicalDevice : devices) {
        VkPhysicalDeviceProperties properties{};
        vkGetPhysicalDeviceProperties(physicalDevice, &properties);
        if (properties.vendorID != kAmdVendorId) {
            continue;
        }

        std::vector<std::string> missing;
        const bool supportsVulkan12 = versionAtLeast(properties.apiVersion, 1U, 2U);
        if (!supportsVulkan12) {
            missing.emplace_back("Vulkan 1.2");
        }
        const std::vector<VkExtensionProperties> extensions =
            enumerateDeviceExtensions(physicalDevice);
        constexpr const char* kPortabilitySubsetExtension =
            "VK_KHR_portability_subset";
        const bool hasPortabilitySubset =
            hasExtension(extensions, kPortabilitySubsetExtension);
        const std::array<const char*, 3> requiredExtensions{
            VK_KHR_ACCELERATION_STRUCTURE_EXTENSION_NAME,
            VK_KHR_DEFERRED_HOST_OPERATIONS_EXTENSION_NAME,
            VK_KHR_RAY_QUERY_EXTENSION_NAME,
        };
        bool hasAllRequiredExtensions = true;
        for (const char* required : requiredExtensions) {
            if (!hasExtension(extensions, required)) {
                missing.emplace_back(required);
                hasAllRequiredExtensions = false;
            }
        }
        const auto [hasComputeQueue, queueFamily] = findComputeQueue(physicalDevice);
        if (!hasComputeQueue) {
            missing.emplace_back("compute queue");
        }

        // Extension feature structures are legal in this query only after their extensions
        // and Vulkan 1.2 have been established for this candidate.
        if (supportsVulkan12 && hasAllRequiredExtensions) {
            VkPhysicalDeviceRayQueryFeaturesKHR rayQueryFeatures{
                VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_RAY_QUERY_FEATURES_KHR};
            VkPhysicalDeviceAccelerationStructureFeaturesKHR accelerationFeatures{
                VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ACCELERATION_STRUCTURE_FEATURES_KHR};
            accelerationFeatures.pNext = &rayQueryFeatures;
            VkPhysicalDeviceBufferDeviceAddressFeatures addressFeatures{
                VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_BUFFER_DEVICE_ADDRESS_FEATURES};
            addressFeatures.pNext = &accelerationFeatures;
            VkPhysicalDeviceFeatures2 features{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2};
            features.pNext = &addressFeatures;
            vkGetPhysicalDeviceFeatures2(physicalDevice, &features);
            if (addressFeatures.bufferDeviceAddress != VK_TRUE) {
                missing.emplace_back("buffer device address feature");
            }
            if (accelerationFeatures.accelerationStructure != VK_TRUE) {
                missing.emplace_back("acceleration structure feature");
            }
            if (rayQueryFeatures.rayQuery != VK_TRUE) {
                missing.emplace_back("ray query feature");
            }
        }

        if (!missing.empty()) {
            std::ostringstream description;
            description << properties.deviceName << " (missing ";
            for (std::size_t index = 0U; index < missing.size(); ++index) {
                if (index != 0U) {
                    description << ", ";
                }
                description << missing[index];
            }
            description << ')';
            rejectedAmdDevices.push_back(description.str());
            continue;
        }

        const int score =
            properties.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU ? 2 : 1;
        if (score <= selectedScore) {
            continue;
        }
        VkPhysicalDeviceAccelerationStructurePropertiesKHR accelerationProperties{
            VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ACCELERATION_STRUCTURE_PROPERTIES_KHR};
        VkPhysicalDeviceProperties2 extendedProperties{
            VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2};
        extendedProperties.pNext = &accelerationProperties;
        vkGetPhysicalDeviceProperties2(physicalDevice, &extendedProperties);

        selected.physicalDevice = physicalDevice;
        selected.properties = properties;
        selected.accelerationProperties = accelerationProperties;
        selected.queueFamily = queueFamily;
        selected.timestampValidBits =
            queueTimestampValidBits(physicalDevice, queueFamily);
        selected.deviceName = properties.deviceName;
        selected.requiresPortabilitySubset = hasPortabilitySubset;
        selected.scratchAlignment = std::max<VkDeviceSize>(
            accelerationProperties.minAccelerationStructureScratchOffsetAlignment, 1U);
        vkGetPhysicalDeviceMemoryProperties(
            physicalDevice, &selected.memoryProperties);
        selectedScore = score;
    }

    if (selected.physicalDevice == VK_NULL_HANDLE) {
        std::ostringstream reason;
        reason << "No AMD Vulkan device satisfies the Ray Query backend requirements.";
        if (!rejectedAmdDevices.empty()) {
            reason << " Rejected devices: ";
            for (std::size_t index = 0U; index < rejectedAmdDevices.size(); ++index) {
                if (index != 0U) {
                    reason << "; ";
                }
                reason << rejectedAmdDevices[index];
            }
            reason << '.';
        }
        throw std::runtime_error(reason.str());
    }
    return selected;
}

template <typename Function>
Function loadDeviceFunction(VkDevice device, const char* name) {
    const auto function = reinterpret_cast<Function>(vkGetDeviceProcAddr(device, name));
    if (function == nullptr) {
        throw std::runtime_error(
            std::string{"Missing Vulkan device function: "} + name + '.');
    }
    return function;
}

struct AccelerationFunctions {
    PFN_vkCreateAccelerationStructureKHR create{nullptr};
    PFN_vkDestroyAccelerationStructureKHR destroy{nullptr};
    PFN_vkGetAccelerationStructureBuildSizesKHR getBuildSizes{nullptr};
    PFN_vkCmdBuildAccelerationStructuresKHR cmdBuild{nullptr};
    PFN_vkGetAccelerationStructureDeviceAddressKHR getDeviceAddress{nullptr};
};

struct LogicalDevice {
    UniqueVulkanHandle<VkDevice> owner;
    VkQueue queue{VK_NULL_HANDLE};
    AccelerationFunctions acceleration;
};

LogicalDevice createLogicalDevice(const PhysicalDeviceSelection& physical) {
    constexpr float kQueuePriority = 1.0F;
    const VkDeviceQueueCreateInfo queueCreateInfo{
        VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO,
        nullptr,
        0U,
        physical.queueFamily,
        1U,
        &kQueuePriority,
    };
    VkPhysicalDeviceRayQueryFeaturesKHR rayQueryFeatures{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_RAY_QUERY_FEATURES_KHR};
    rayQueryFeatures.rayQuery = VK_TRUE;
    VkPhysicalDeviceAccelerationStructureFeaturesKHR accelerationFeatures{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ACCELERATION_STRUCTURE_FEATURES_KHR};
    accelerationFeatures.pNext = &rayQueryFeatures;
    accelerationFeatures.accelerationStructure = VK_TRUE;
    VkPhysicalDeviceBufferDeviceAddressFeatures addressFeatures{
        VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_BUFFER_DEVICE_ADDRESS_FEATURES};
    addressFeatures.pNext = &accelerationFeatures;
    addressFeatures.bufferDeviceAddress = VK_TRUE;
    constexpr const char* kPortabilitySubsetExtension = "VK_KHR_portability_subset";
    const std::array<const char*, 4> extensions{
        VK_KHR_ACCELERATION_STRUCTURE_EXTENSION_NAME,
        VK_KHR_DEFERRED_HOST_OPERATIONS_EXTENSION_NAME,
        VK_KHR_RAY_QUERY_EXTENSION_NAME,
        kPortabilitySubsetExtension,
    };
    const std::uint32_t extensionCount =
        physical.requiresPortabilitySubset ? 4U : 3U;
    const VkDeviceCreateInfo createInfo{
        VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO,
        &addressFeatures,
        0U,
        1U,
        &queueCreateInfo,
        0U,
        nullptr,
        extensionCount,
        extensions.data(),
        nullptr,
    };

    VkDevice device = VK_NULL_HANDLE;
    requireVulkanSuccess(
        vkCreateDevice(physical.physicalDevice, &createInfo, nullptr, &device),
        "vkCreateDevice");
    LogicalDevice result;
    result.owner.reset(device, [](VkDevice value) { vkDestroyDevice(value, nullptr); });
    vkGetDeviceQueue(device, physical.queueFamily, 0U, &result.queue);
    result.acceleration.create =
        loadDeviceFunction<PFN_vkCreateAccelerationStructureKHR>(
            device, "vkCreateAccelerationStructureKHR");
    result.acceleration.destroy =
        loadDeviceFunction<PFN_vkDestroyAccelerationStructureKHR>(
            device, "vkDestroyAccelerationStructureKHR");
    result.acceleration.getBuildSizes =
        loadDeviceFunction<PFN_vkGetAccelerationStructureBuildSizesKHR>(
            device, "vkGetAccelerationStructureBuildSizesKHR");
    result.acceleration.cmdBuild =
        loadDeviceFunction<PFN_vkCmdBuildAccelerationStructuresKHR>(
            device, "vkCmdBuildAccelerationStructuresKHR");
    result.acceleration.getDeviceAddress =
        loadDeviceFunction<PFN_vkGetAccelerationStructureDeviceAddressKHR>(
            device, "vkGetAccelerationStructureDeviceAddressKHR");
    return result;
}

std::uint32_t findMemoryType(
    const VkPhysicalDeviceMemoryProperties& properties,
    std::uint32_t allowedTypes,
    VkMemoryPropertyFlags requiredProperties) {
    for (std::uint32_t index = 0U; index < properties.memoryTypeCount; ++index) {
        const bool allowed = (allowedTypes & (1U << index)) != 0U;
        const bool satisfies =
            (properties.memoryTypes[index].propertyFlags & requiredProperties) ==
            requiredProperties;
        if (allowed && satisfies) {
            return index;
        }
    }
    throw std::runtime_error(
        "No Vulkan memory type satisfies the requested properties.");
}

}  // namespace

class VulkanRuntime::Implementation {
public:
    Implementation(const PackedSceneData& scene, VulkanRayQueryOptions options);
    ~Implementation();

    [[nodiscard]] const std::string& deviceName() const noexcept;
    [[nodiscard]] const VulkanRayQuerySetupTimings& setupTimings() const noexcept;
    [[nodiscard]] VulkanValidationReport validationReport() const;
    [[nodiscard]] std::size_t triangleCount() const noexcept;
    [[nodiscard]] VkDevice device() const noexcept;
    [[nodiscard]] const VkPhysicalDeviceProperties& physicalProperties() const noexcept;
    [[nodiscard]] const VkPhysicalDeviceMemoryProperties& memoryProperties() const noexcept;
    [[nodiscard]] std::uint32_t timestampValidBits() const noexcept;
    [[nodiscard]] VkAccelerationStructureKHR topLevelAccelerationStructure() const noexcept;
    [[nodiscard]] const VulkanBuffer& vertexBuffer() const noexcept;
    [[nodiscard]] const VulkanBuffer& indexBuffer() const noexcept;
    [[nodiscard]] const VulkanBuffer& triangleMaterialIdBuffer() const noexcept;
    [[nodiscard]] const VulkanBuffer& materialBuffer() const noexcept;
    [[nodiscard]] const VulkanBuffer& textureDescriptorBuffer() const noexcept;
    [[nodiscard]] const VulkanBuffer& textureMipBuffer() const noexcept;
    [[nodiscard]] const VulkanBuffer& textureTexelBuffer() const noexcept;
    [[nodiscard]] const VulkanBuffer& areaLightBuffer() const noexcept;
    [[nodiscard]] const VulkanBuffer& areaLightTriangleBuffer() const noexcept;
    [[nodiscard]] const VulkanBuffer& environmentRowBuffer() const noexcept;
    [[nodiscard]] const VulkanBuffer& environmentTexelBuffer() const noexcept;
    [[nodiscard]] std::mutex& operationMutex() noexcept;
    [[nodiscard]] VkCommandBuffer beginCommands();
    void submitAndWait(const char* description);

private:
    struct State;

    void createCommandResources();
    void createGeometry(const PackedSceneData& scene);

    std::unique_ptr<State> state_;
};

double elapsedMilliseconds(VulkanClock::time_point begin) noexcept {
    return std::chrono::duration<double, std::milli>(VulkanClock::now() - begin).count();
}

void requireVulkanSuccess(VkResult result, const char* operation) {
    if (result != VK_SUCCESS) {
        throw std::runtime_error(
            std::string{operation} + " failed with Vulkan result " +
            std::to_string(static_cast<int>(result)) + '.');
    }
}

std::vector<std::uint32_t> readSpirvFile(const std::filesystem::path& path) {
    std::ifstream stream{path, std::ios::binary | std::ios::ate};
    if (!stream) {
        throw std::runtime_error("Unable to open SPIR-V file: " + path.string());
    }
    const std::streamoff end = stream.tellg();
    if (end <= 0 ||
        static_cast<std::uint64_t>(end) % sizeof(std::uint32_t) != 0U ||
        static_cast<std::uint64_t>(end) > std::numeric_limits<std::size_t>::max() ||
        static_cast<std::uint64_t>(end) >
            static_cast<std::uint64_t>(std::numeric_limits<std::streamsize>::max())) {
        throw std::runtime_error(
            "The SPIR-V file has an invalid size: " + path.string());
    }
    const std::size_t byteCount = static_cast<std::size_t>(end);
    std::vector<std::uint32_t> words(byteCount / sizeof(std::uint32_t));
    stream.seekg(0, std::ios::beg);
    stream.read(
        reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(byteCount));
    if (!stream) {
        throw std::runtime_error(
            "Unable to read the complete SPIR-V file: " + path.string());
    }
    constexpr std::uint32_t kSpirvMagic = 0x07230203U;
    constexpr std::uint32_t kSpirv10 = 0x00010000U;
    constexpr std::uint32_t kSpirv15 = 0x00010500U;
    if (words.size() < 5U || words[0] != kSpirvMagic || words[1] < kSpirv10 ||
        words[1] > kSpirv15 || words[3] == 0U || words[4] != 0U) {
        throw std::runtime_error(
            "The file does not contain SPIR-V bytecode: " + path.string());
    }
    return words;
}

VkDeviceSize checkedVulkanByteSize(
    std::size_t count, std::size_t stride, const char* description) {
    if (count == 0U || stride == 0U ||
        count > std::numeric_limits<std::size_t>::max() / stride) {
        throw std::overflow_error(
            std::string{description} + " byte size overflow.");
    }
    const std::size_t bytes = count * stride;
    if (bytes > std::numeric_limits<VkDeviceSize>::max()) {
        throw std::overflow_error(
            std::string{description} + " exceeds VkDeviceSize.");
    }
    return static_cast<VkDeviceSize>(bytes);
}

class VulkanBuffer::Implementation {
public:
    Implementation(
        VkDevice device,
        const VkPhysicalDeviceMemoryProperties& memoryProperties,
        VkDeviceSize size,
        VkBufferUsageFlags usage,
        VkMemoryPropertyFlags requiredMemoryProperties,
        bool needsDeviceAddress)
        : device_(device), size_(size) {
        if (size == 0U) {
            throw std::invalid_argument("A zero-sized Vulkan buffer was requested.");
        }
        const VkBufferCreateInfo bufferCreateInfo{
            VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO,
            nullptr,
            0U,
            size,
            usage,
            VK_SHARING_MODE_EXCLUSIVE,
            0U,
            nullptr,
        };
        VkBuffer buffer = VK_NULL_HANDLE;
        requireVulkanSuccess(
            vkCreateBuffer(device, &bufferCreateInfo, nullptr, &buffer),
            "vkCreateBuffer");
        buffer_.reset(
            buffer,
            [device](VkBuffer value) { vkDestroyBuffer(device, value, nullptr); });

        VkMemoryRequirements requirements{};
        vkGetBufferMemoryRequirements(device, buffer, &requirements);
        const VkMemoryAllocateFlagsInfo flagsInfo{
            VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_FLAGS_INFO,
            nullptr,
            needsDeviceAddress ? VK_MEMORY_ALLOCATE_DEVICE_ADDRESS_BIT : 0U,
            0U,
        };
        const VkMemoryAllocateInfo allocationInfo{
            VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO,
            needsDeviceAddress ? &flagsInfo : nullptr,
            requirements.size,
            findMemoryType(
                memoryProperties,
                requirements.memoryTypeBits,
                requiredMemoryProperties),
        };
        VkDeviceMemory memory = VK_NULL_HANDLE;
        requireVulkanSuccess(
            vkAllocateMemory(device, &allocationInfo, nullptr, &memory),
            "vkAllocateMemory");
        memory_.reset(
            memory,
            [device](VkDeviceMemory value) { vkFreeMemory(device, value, nullptr); });
        requireVulkanSuccess(
            vkBindBufferMemory(device, buffer, memory, 0U), "vkBindBufferMemory");
    }

    void write(const void* source, VkDeviceSize byteCount) const {
        if (source == nullptr || byteCount == 0U || byteCount > size_ ||
            byteCount > std::numeric_limits<std::size_t>::max()) {
            throw std::invalid_argument("Invalid Vulkan buffer write.");
        }
        void* mapped = nullptr;
        requireVulkanSuccess(
            vkMapMemory(device_, memory_.get(), 0U, byteCount, 0U, &mapped),
            "vkMapMemory(write)");
        std::memcpy(mapped, source, static_cast<std::size_t>(byteCount));
        vkUnmapMemory(device_, memory_.get());
    }

    void read(void* destination, VkDeviceSize byteCount) const {
        if (destination == nullptr || byteCount == 0U || byteCount > size_ ||
            byteCount > std::numeric_limits<std::size_t>::max()) {
            throw std::invalid_argument("Invalid Vulkan buffer read.");
        }
        void* mapped = nullptr;
        requireVulkanSuccess(
            vkMapMemory(device_, memory_.get(), 0U, byteCount, 0U, &mapped),
            "vkMapMemory(read)");
        std::memcpy(destination, mapped, static_cast<std::size_t>(byteCount));
        vkUnmapMemory(device_, memory_.get());
    }

    [[nodiscard]] VkBuffer get() const noexcept {
        return buffer_.get();
    }

    [[nodiscard]] VkDeviceSize size() const noexcept {
        return size_;
    }

    [[nodiscard]] VkDeviceAddress deviceAddress() const {
        const VkBufferDeviceAddressInfo info{
            VK_STRUCTURE_TYPE_BUFFER_DEVICE_ADDRESS_INFO,
            nullptr,
            buffer_.get(),
        };
        const VkDeviceAddress result = vkGetBufferDeviceAddress(device_, &info);
        if (result == 0U) {
            throw std::runtime_error(
                "vkGetBufferDeviceAddress returned a null address.");
        }
        return result;
    }

private:
    VkDevice device_{VK_NULL_HANDLE};
    VkDeviceSize size_{0U};
    UniqueVulkanHandle<VkDeviceMemory> memory_;
    UniqueVulkanHandle<VkBuffer> buffer_;
};

VulkanBuffer::VulkanBuffer(
    VkDevice device,
    const VkPhysicalDeviceMemoryProperties& memoryProperties,
    VkDeviceSize size,
    VkBufferUsageFlags usage,
    VkMemoryPropertyFlags requiredMemoryProperties,
    bool needsDeviceAddress)
    : implementation_(std::make_unique<Implementation>(
          device,
          memoryProperties,
          size,
          usage,
          requiredMemoryProperties,
          needsDeviceAddress)) {}

VulkanBuffer::~VulkanBuffer() = default;

void VulkanBuffer::write(const void* source, VkDeviceSize byteCount) const {
    implementation_->write(source, byteCount);
}

void VulkanBuffer::read(void* destination, VkDeviceSize byteCount) const {
    implementation_->read(destination, byteCount);
}

VkBuffer VulkanBuffer::get() const noexcept {
    return implementation_->get();
}

VkDeviceSize VulkanBuffer::size() const noexcept {
    return implementation_->size();
}

VkDeviceAddress VulkanBuffer::deviceAddress() const {
    return implementation_->deviceAddress();
}

namespace {

class AccelerationStructure {
public:
    AccelerationStructure(
        VkDevice device,
        const AccelerationFunctions& functions,
        const VkPhysicalDeviceMemoryProperties& memoryProperties,
        VkDeviceSize size,
        VkAccelerationStructureTypeKHR type)
        : device_(device),
          functions_(functions),
          storage_(
              device,
              memoryProperties,
              size,
              VK_BUFFER_USAGE_ACCELERATION_STRUCTURE_STORAGE_BIT_KHR |
                  VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT,
              VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
              true) {
        const VkAccelerationStructureCreateInfoKHR info{
            VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_CREATE_INFO_KHR,
            nullptr,
            0U,
            storage_.get(),
            0U,
            size,
            type,
            0U,
        };
        requireVulkanSuccess(
            functions_.create(device_, &info, nullptr, &handle_),
            "vkCreateAccelerationStructureKHR");
    }

    AccelerationStructure(const AccelerationStructure&) = delete;
    AccelerationStructure& operator=(const AccelerationStructure&) = delete;

    ~AccelerationStructure() {
        if (handle_ != VK_NULL_HANDLE) {
            functions_.destroy(device_, handle_, nullptr);
        }
    }

    [[nodiscard]] VkAccelerationStructureKHR get() const noexcept {
        return handle_;
    }

    [[nodiscard]] VkDeviceAddress deviceAddress() const {
        const VkAccelerationStructureDeviceAddressInfoKHR info{
            VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_DEVICE_ADDRESS_INFO_KHR,
            nullptr,
            handle_,
        };
        const VkDeviceAddress result =
            functions_.getDeviceAddress(device_, &info);
        if (result == 0U) {
            throw std::runtime_error(
                "vkGetAccelerationStructureDeviceAddressKHR returned a null address.");
        }
        return result;
    }

private:
    VkDevice device_{VK_NULL_HANDLE};
    const AccelerationFunctions& functions_;
    VulkanBuffer storage_;
    VkAccelerationStructureKHR handle_{VK_NULL_HANDLE};
};

VkDeviceSize checkedScratchSize(VkDeviceSize requiredSize, VkDeviceSize alignment) {
    requiredSize = std::max<VkDeviceSize>(requiredSize, 1U);
    if (alignment == 0U ||
        requiredSize > std::numeric_limits<VkDeviceSize>::max() - (alignment - 1U)) {
        throw std::overflow_error("Acceleration-structure scratch size overflow.");
    }
    return requiredSize + alignment - 1U;
}

VkDeviceAddress alignAddress(VkDeviceAddress address, VkDeviceSize alignment) {
    const VkDeviceAddress remainder = address % alignment;
    if (remainder == 0U) {
        return address;
    }
    const VkDeviceAddress adjustment = alignment - remainder;
    if (address > std::numeric_limits<VkDeviceAddress>::max() - adjustment) {
        throw std::overflow_error("Acceleration-structure scratch address overflow.");
    }
    return address + adjustment;
}

bool hasReferencedCutout(const PackedSceneData& scene) noexcept {
    return std::any_of(
        scene.triangleMaterialIds.begin(),
        scene.triangleMaterialIds.end(),
        [&scene](std::uint32_t materialId) {
            return (scene.materials[materialId].flagsAndReserved[0] &
                    kPackedMaterialCutout) != 0U;
        });
}

template <typename Element>
VkDeviceSize packedArrayBufferByteSize(
    const std::vector<Element>& elements, const char* description) {
    return elements.empty()
               ? static_cast<VkDeviceSize>(sizeof(Element))
               : checkedVulkanByteSize(elements.size(), sizeof(Element), description);
}

void validatePackedSceneForGpu(
    const PackedSceneData& scene, const PhysicalDeviceSelection& physical) {
    const std::size_t triangleCount = scene.triangleCount();
    if (triangleCount > physical.accelerationProperties.maxPrimitiveCount ||
        physical.accelerationProperties.maxGeometryCount < 1U ||
        physical.accelerationProperties.maxInstanceCount < 1U) {
        throw std::invalid_argument(
            "The packed scene exceeds this device's acceleration-structure limits.");
    }
    for (std::uint32_t materialId : scene.triangleMaterialIds) {
        const PackedMaterial& material = scene.materials[materialId];
        if ((material.flagsAndReserved[0] & kPackedMaterialCutout) != 0U &&
            ((material.flagsAndReserved[0] &
              kPackedMaterialHasDiffuseTexture) == 0U ||
             material.textureIds[0] == kInvalidSceneId)) {
            throw std::invalid_argument(
                "Vulkan alpha-cutout materials require a packed diffuse texture.");
        }
    }
    const VkDeviceSize vertexBytes = checkedVulkanByteSize(
        scene.vertices.size(), sizeof(PackedVertex), "Packed vertex buffer");
    const VkDeviceSize indexBytes = checkedVulkanByteSize(
        scene.triangleIndices.size(),
        sizeof(std::uint32_t),
        "Packed index buffer");
    const VkDeviceSize triangleMaterialIdBytes = checkedVulkanByteSize(
        scene.triangleMaterialIds.size(),
        sizeof(std::uint32_t),
        "Packed triangle material ID buffer");
    const VkDeviceSize materialBytes = checkedVulkanByteSize(
        scene.materials.size(), sizeof(PackedMaterial), "Packed material buffer");

    const VkDeviceSize textureDescriptorBytes = packedArrayBufferByteSize(
        scene.textures, "Packed texture descriptor buffer");
    const VkDeviceSize textureMipBytes = packedArrayBufferByteSize(
        scene.textureMipLevels, "Packed texture mip buffer");
    const VkDeviceSize textureTexelBytes = packedArrayBufferByteSize(
        scene.textureTexelsRgba8, "Packed texture texel buffer");
    const VkDeviceSize areaLightBytes = packedArrayBufferByteSize(
        scene.areaLights, "Packed area light buffer");
    const VkDeviceSize areaLightTriangleBytes = packedArrayBufferByteSize(
        scene.areaLightTriangles, "Packed area light triangle buffer");
    const VkDeviceSize environmentRowBytes = packedArrayBufferByteSize(
        scene.environmentRows, "Packed environment row buffer");
    const VkDeviceSize environmentTexelBytes = packedArrayBufferByteSize(
        scene.environmentTexels, "Packed environment texel buffer");
    const VkDeviceSize maximumStorageRange =
        physical.properties.limits.maxStorageBufferRange;
    const auto requireStorageBufferRange =
        [maximumStorageRange](VkDeviceSize bytes, const char* description) {
            if (bytes > maximumStorageRange) {
                throw std::invalid_argument(
                    std::string{description} + " exceeds maxStorageBufferRange.");
            }
        };
    requireStorageBufferRange(vertexBytes, "Packed vertex buffer");
    requireStorageBufferRange(indexBytes, "Packed index buffer");
    requireStorageBufferRange(
        triangleMaterialIdBytes, "Packed triangle material ID buffer");
    requireStorageBufferRange(materialBytes, "Packed material buffer");
    requireStorageBufferRange(
        textureDescriptorBytes, "Packed texture descriptor buffer");
    requireStorageBufferRange(textureMipBytes, "Packed texture mip buffer");
    requireStorageBufferRange(textureTexelBytes, "Packed texture texel buffer");
    requireStorageBufferRange(areaLightBytes, "Packed area light buffer");
    requireStorageBufferRange(
        areaLightTriangleBytes, "Packed area light triangle buffer");
    requireStorageBufferRange(environmentRowBytes, "Packed environment row buffer");
    requireStorageBufferRange(
        environmentTexelBytes, "Packed environment texel buffer");
}

}  // namespace

struct VulkanRuntime::Implementation::State {
    State(const PackedSceneData& scene, VulkanRayQueryOptions runtimeOptions)
        : validation(runtimeOptions.requestValidation),
          instance(runtimeOptions.requestValidation, validation),
          physical(selectAmdDevice(instance.get())),
          logical(createLogicalDevice(physical)),
          options(runtimeOptions),
          triangleCount(scene.triangleCount()),
          hasAlphaCutout(hasReferencedCutout(scene)) {}

    ValidationState validation;
    VulkanInstance instance;
    PhysicalDeviceSelection physical;
    LogicalDevice logical;
    VulkanRayQueryOptions options;
    VulkanRayQuerySetupTimings timings;
    std::size_t triangleCount{0U};
    bool hasAlphaCutout{false};
    std::mutex operationMutex;

    std::unique_ptr<VulkanBuffer> vertexBuffer;
    std::unique_ptr<VulkanBuffer> indexBuffer;
    std::unique_ptr<VulkanBuffer> triangleMaterialIdBuffer;
    std::unique_ptr<VulkanBuffer> materialBuffer;
    std::unique_ptr<VulkanBuffer> textureDescriptorBuffer;
    std::unique_ptr<VulkanBuffer> textureMipBuffer;
    std::unique_ptr<VulkanBuffer> textureTexelBuffer;
    std::unique_ptr<VulkanBuffer> areaLightBuffer;
    std::unique_ptr<VulkanBuffer> areaLightTriangleBuffer;
    std::unique_ptr<VulkanBuffer> environmentRowBuffer;
    std::unique_ptr<VulkanBuffer> environmentTexelBuffer;
    std::unique_ptr<VulkanBuffer> instanceBuffer;
    std::unique_ptr<VulkanBuffer> scratchBuffer;
    std::unique_ptr<AccelerationStructure> bottomLevel;
    std::unique_ptr<AccelerationStructure> topLevel;

    UniqueVulkanHandle<VkCommandPool> commandPool;
    VkCommandBuffer commandBuffer{VK_NULL_HANDLE};
    UniqueVulkanHandle<VkFence> fence;
};

VulkanRuntime::Implementation::Implementation(
    const PackedSceneData& scene, VulkanRayQueryOptions options)
    : state_(std::make_unique<State>(scene, options)) {
    if (state_->options.fenceTimeoutNanoseconds == 0U) {
        throw std::invalid_argument("The Vulkan fence timeout must be non-zero.");
    }
    validatePackedSceneForGpu(scene, state_->physical);
    createCommandResources();
    createGeometry(scene);
}

VulkanRuntime::Implementation::~Implementation() {
    if (state_->logical.owner.get() != VK_NULL_HANDLE) {
        // Public operations are synchronous. This is a final defensive quiescence before
        // dependent resources begin their reverse-order destruction.
        (void)vkDeviceWaitIdle(state_->logical.owner.get());
    }
}

const std::string& VulkanRuntime::Implementation::deviceName() const noexcept {
    return state_->physical.deviceName;
}

const VulkanRayQuerySetupTimings&
VulkanRuntime::Implementation::setupTimings() const noexcept {
    return state_->timings;
}

VulkanValidationReport VulkanRuntime::Implementation::validationReport() const {
    return state_->validation.snapshot();
}

std::size_t VulkanRuntime::Implementation::triangleCount() const noexcept {
    return state_->triangleCount;
}

VkDevice VulkanRuntime::Implementation::device() const noexcept {
    return state_->logical.owner.get();
}

const VkPhysicalDeviceProperties&
VulkanRuntime::Implementation::physicalProperties() const noexcept {
    return state_->physical.properties;
}

const VkPhysicalDeviceMemoryProperties&
VulkanRuntime::Implementation::memoryProperties() const noexcept {
    return state_->physical.memoryProperties;
}

std::uint32_t VulkanRuntime::Implementation::timestampValidBits() const noexcept {
    return state_->physical.timestampValidBits;
}

VkAccelerationStructureKHR
VulkanRuntime::Implementation::topLevelAccelerationStructure() const noexcept {
    return state_->topLevel->get();
}

const VulkanBuffer& VulkanRuntime::Implementation::vertexBuffer() const noexcept {
    return *state_->vertexBuffer;
}

const VulkanBuffer& VulkanRuntime::Implementation::indexBuffer() const noexcept {
    return *state_->indexBuffer;
}

const VulkanBuffer&
VulkanRuntime::Implementation::triangleMaterialIdBuffer() const noexcept {
    return *state_->triangleMaterialIdBuffer;
}

const VulkanBuffer& VulkanRuntime::Implementation::materialBuffer() const noexcept {
    return *state_->materialBuffer;
}

const VulkanBuffer&
VulkanRuntime::Implementation::textureDescriptorBuffer() const noexcept {
    return *state_->textureDescriptorBuffer;
}

const VulkanBuffer& VulkanRuntime::Implementation::textureMipBuffer() const noexcept {
    return *state_->textureMipBuffer;
}

const VulkanBuffer& VulkanRuntime::Implementation::textureTexelBuffer() const noexcept {
    return *state_->textureTexelBuffer;
}

const VulkanBuffer& VulkanRuntime::Implementation::areaLightBuffer() const noexcept {
    return *state_->areaLightBuffer;
}

const VulkanBuffer&
VulkanRuntime::Implementation::areaLightTriangleBuffer() const noexcept {
    return *state_->areaLightTriangleBuffer;
}

const VulkanBuffer&
VulkanRuntime::Implementation::environmentRowBuffer() const noexcept {
    return *state_->environmentRowBuffer;
}

const VulkanBuffer&
VulkanRuntime::Implementation::environmentTexelBuffer() const noexcept {
    return *state_->environmentTexelBuffer;
}

std::mutex& VulkanRuntime::Implementation::operationMutex() noexcept {
    return state_->operationMutex;
}

VkCommandBuffer VulkanRuntime::Implementation::beginCommands() {
    const VkDevice logicalDevice = device();
    requireVulkanSuccess(
        vkResetCommandPool(logicalDevice, state_->commandPool.get(), 0U),
        "vkResetCommandPool");
    const VkCommandBufferBeginInfo beginInfo{
        VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
        nullptr,
        VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT,
        nullptr,
    };
    requireVulkanSuccess(
        vkBeginCommandBuffer(state_->commandBuffer, &beginInfo),
        "vkBeginCommandBuffer");
    return state_->commandBuffer;
}

void VulkanRuntime::Implementation::submitAndWait(const char* description) {
    const VkDevice logicalDevice = device();
    requireVulkanSuccess(
        vkEndCommandBuffer(state_->commandBuffer), "vkEndCommandBuffer");
    const VkFence fence = state_->fence.get();
    requireVulkanSuccess(
        vkResetFences(logicalDevice, 1U, &fence), "vkResetFences");
    const VkSubmitInfo submitInfo{
        VK_STRUCTURE_TYPE_SUBMIT_INFO,
        nullptr,
        0U,
        nullptr,
        nullptr,
        1U,
        &state_->commandBuffer,
        0U,
        nullptr,
    };
    requireVulkanSuccess(
        vkQueueSubmit(
            state_->logical.queue, 1U, &submitInfo, state_->fence.get()),
        "vkQueueSubmit");
    const VkResult waitResult = vkWaitForFences(
        logicalDevice,
        1U,
        &fence,
        VK_TRUE,
        state_->options.fenceTimeoutNanoseconds);
    if (waitResult != VK_SUCCESS) {
        // Once submitted, no dependent resource may unwind while the queue can still use it.
        const VkResult queueIdle = vkQueueWaitIdle(state_->logical.queue);
        if (queueIdle != VK_SUCCESS && queueIdle != VK_ERROR_DEVICE_LOST) {
            const VkResult deviceIdle = vkDeviceWaitIdle(logicalDevice);
            if (deviceIdle != VK_SUCCESS && deviceIdle != VK_ERROR_DEVICE_LOST) {
                throw std::runtime_error(
                    std::string{description} +
                    " could not quiesce the Vulkan device after a fence failure.");
            }
        }
        if (waitResult == VK_TIMEOUT) {
            throw std::runtime_error(
                std::string{description} + " exceeded its Vulkan fence timeout.");
        }
        requireVulkanSuccess(waitResult, "vkWaitForFences");
    }
}

void VulkanRuntime::Implementation::createCommandResources() {
    const VkDevice logicalDevice = device();
    const VkCommandPoolCreateInfo poolInfo{
        VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,
        nullptr,
        VK_COMMAND_POOL_CREATE_TRANSIENT_BIT |
            VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT,
        state_->physical.queueFamily,
    };
    VkCommandPool commandPool = VK_NULL_HANDLE;
    requireVulkanSuccess(
        vkCreateCommandPool(logicalDevice, &poolInfo, nullptr, &commandPool),
        "vkCreateCommandPool");
    state_->commandPool.reset(
        commandPool,
        [logicalDevice](VkCommandPool value) {
            vkDestroyCommandPool(logicalDevice, value, nullptr);
        });
    const VkCommandBufferAllocateInfo allocateInfo{
        VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
        nullptr,
        commandPool,
        VK_COMMAND_BUFFER_LEVEL_PRIMARY,
        1U,
    };
    requireVulkanSuccess(
        vkAllocateCommandBuffers(
            logicalDevice, &allocateInfo, &state_->commandBuffer),
        "vkAllocateCommandBuffers");

    const VkFenceCreateInfo fenceInfo{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    VkFence fence = VK_NULL_HANDLE;
    requireVulkanSuccess(
        vkCreateFence(logicalDevice, &fenceInfo, nullptr, &fence),
        "vkCreateFence");
    state_->fence.reset(
        fence,
        [logicalDevice](VkFence value) {
            vkDestroyFence(logicalDevice, value, nullptr);
        });
}

void VulkanRuntime::Implementation::createGeometry(const PackedSceneData& scene) {
    const VulkanClock::time_point uploadBegin = VulkanClock::now();
    const VkDevice logicalDevice = device();
    const VkDeviceSize vertexBytes = checkedVulkanByteSize(
        scene.vertices.size(), sizeof(PackedVertex), "Packed vertex buffer");
    const VkDeviceSize indexBytes = checkedVulkanByteSize(
        scene.triangleIndices.size(),
        sizeof(std::uint32_t),
        "Packed index buffer");
    const VkDeviceSize triangleMaterialIdBytes = checkedVulkanByteSize(
        scene.triangleMaterialIds.size(),
        sizeof(std::uint32_t),
        "Packed triangle material ID buffer");
    const VkDeviceSize materialBytes = checkedVulkanByteSize(
        scene.materials.size(), sizeof(PackedMaterial), "Packed material buffer");
    const VkDeviceSize textureDescriptorBytes = packedArrayBufferByteSize(
        scene.textures, "Packed texture descriptor buffer");
    const VkDeviceSize textureMipBytes = packedArrayBufferByteSize(
        scene.textureMipLevels, "Packed texture mip buffer");
    const VkDeviceSize textureTexelBytes = packedArrayBufferByteSize(
        scene.textureTexelsRgba8, "Packed texture texel buffer");
    const VkDeviceSize areaLightBytes = packedArrayBufferByteSize(
        scene.areaLights, "Packed area light buffer");
    const VkDeviceSize areaLightTriangleBytes = packedArrayBufferByteSize(
        scene.areaLightTriangles, "Packed area light triangle buffer");
    const VkDeviceSize environmentRowBytes = packedArrayBufferByteSize(
        scene.environmentRows, "Packed environment row buffer");
    const VkDeviceSize environmentTexelBytes = packedArrayBufferByteSize(
        scene.environmentTexels, "Packed environment texel buffer");
    constexpr VkBufferUsageFlags kGeometryUsage =
        VK_BUFFER_USAGE_TRANSFER_DST_BIT |
        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT |
        VK_BUFFER_USAGE_ACCELERATION_STRUCTURE_BUILD_INPUT_READ_ONLY_BIT_KHR |
        VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT;
    constexpr VkBufferUsageFlags kShadingDataUsage =
        VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_STORAGE_BUFFER_BIT;
    state_->vertexBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        vertexBytes,
        kGeometryUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        true);
    state_->indexBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        indexBytes,
        kGeometryUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        true);
    state_->triangleMaterialIdBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        triangleMaterialIdBytes,
        kShadingDataUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        false);
    state_->materialBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        materialBytes,
        kShadingDataUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        false);
    state_->textureDescriptorBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        textureDescriptorBytes,
        kShadingDataUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        false);
    state_->textureMipBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        textureMipBytes,
        kShadingDataUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        false);
    state_->textureTexelBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        textureTexelBytes,
        kShadingDataUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        false);
    state_->areaLightBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        areaLightBytes,
        kShadingDataUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        false);
    state_->areaLightTriangleBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        areaLightTriangleBytes,
        kShadingDataUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        false);
    state_->environmentRowBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        environmentRowBytes,
        kShadingDataUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        false);
    state_->environmentTexelBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        environmentTexelBytes,
        kShadingDataUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        false);

    VkAccelerationStructureGeometryKHR bottomGeometry{
        VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_GEOMETRY_KHR};
    bottomGeometry.geometryType = VK_GEOMETRY_TYPE_TRIANGLES_KHR;
    bottomGeometry.flags =
        state_->hasAlphaCutout
            ? VK_GEOMETRY_NO_DUPLICATE_ANY_HIT_INVOCATION_BIT_KHR
            : VK_GEOMETRY_OPAQUE_BIT_KHR;
    bottomGeometry.geometry.triangles.sType =
        VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_GEOMETRY_TRIANGLES_DATA_KHR;
    bottomGeometry.geometry.triangles.vertexFormat = VK_FORMAT_R32G32B32_SFLOAT;
    bottomGeometry.geometry.triangles.vertexData.deviceAddress =
        state_->vertexBuffer->deviceAddress();
    if (bottomGeometry.geometry.triangles.vertexData.deviceAddress % sizeof(float) != 0U) {
        throw std::runtime_error("The Vulkan vertex buffer address is misaligned.");
    }
    bottomGeometry.geometry.triangles.vertexStride = sizeof(PackedVertex);
    bottomGeometry.geometry.triangles.maxVertex =
        static_cast<std::uint32_t>(scene.vertices.size() - 1U);
    bottomGeometry.geometry.triangles.indexType = VK_INDEX_TYPE_UINT32;
    bottomGeometry.geometry.triangles.indexData.deviceAddress =
        state_->indexBuffer->deviceAddress();
    if (bottomGeometry.geometry.triangles.indexData.deviceAddress %
            sizeof(std::uint32_t) !=
        0U) {
        throw std::runtime_error("The Vulkan index buffer address is misaligned.");
    }

    VkAccelerationStructureBuildGeometryInfoKHR bottomBuildInfo{
        VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_BUILD_GEOMETRY_INFO_KHR};
    bottomBuildInfo.type = VK_ACCELERATION_STRUCTURE_TYPE_BOTTOM_LEVEL_KHR;
    bottomBuildInfo.flags = VK_BUILD_ACCELERATION_STRUCTURE_PREFER_FAST_TRACE_BIT_KHR;
    bottomBuildInfo.mode = VK_BUILD_ACCELERATION_STRUCTURE_MODE_BUILD_KHR;
    bottomBuildInfo.geometryCount = 1U;
    bottomBuildInfo.pGeometries = &bottomGeometry;
    const std::uint32_t primitiveCount =
        static_cast<std::uint32_t>(state_->triangleCount);
    VkAccelerationStructureBuildSizesInfoKHR bottomSizes{
        VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_BUILD_SIZES_INFO_KHR};
    state_->logical.acceleration.getBuildSizes(
        logicalDevice,
        VK_ACCELERATION_STRUCTURE_BUILD_TYPE_DEVICE_KHR,
        &bottomBuildInfo,
        &primitiveCount,
        &bottomSizes);
    state_->bottomLevel = std::make_unique<AccelerationStructure>(
        logicalDevice,
        state_->logical.acceleration,
        state_->physical.memoryProperties,
        bottomSizes.accelerationStructureSize,
        VK_ACCELERATION_STRUCTURE_TYPE_BOTTOM_LEVEL_KHR);
    bottomBuildInfo.dstAccelerationStructure = state_->bottomLevel->get();

    VkAccelerationStructureInstanceKHR instanceRecord{};
    instanceRecord.transform.matrix[0][0] = 1.0F;
    instanceRecord.transform.matrix[1][1] = 1.0F;
    instanceRecord.transform.matrix[2][2] = 1.0F;
    instanceRecord.instanceCustomIndex = 0U;
    instanceRecord.mask = kInstanceMask;
    instanceRecord.instanceShaderBindingTableRecordOffset = 0U;
    instanceRecord.flags =
        VK_GEOMETRY_INSTANCE_TRIANGLE_FACING_CULL_DISABLE_BIT_KHR;
    instanceRecord.accelerationStructureReference =
        state_->bottomLevel->deviceAddress();
    const VkDeviceSize instanceBytes = sizeof(instanceRecord);
    state_->instanceBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        instanceBytes,
        kGeometryUsage,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        true);
    if (state_->instanceBuffer->deviceAddress() % 16U != 0U) {
        throw std::runtime_error("The Vulkan instance buffer address is misaligned.");
    }

    constexpr VkMemoryPropertyFlags kHostUploadMemory =
        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
        VK_MEMORY_PROPERTY_HOST_COHERENT_BIT;
    const PackedTexture dummyTexture{};
    const PackedTextureMip dummyMip{};
    const std::uint32_t dummyTexel = 0U;
    const PackedAreaLight dummyAreaLight{};
    const PackedAreaLightTriangle dummyAreaLightTriangle{};
    const PackedEnvironmentRow dummyEnvironmentRow{};
    const PackedEnvironmentTexel dummyEnvironmentTexel{};
    const void* textureDescriptorSource =
        scene.textures.empty()
            ? static_cast<const void*>(&dummyTexture)
            : static_cast<const void*>(scene.textures.data());
    const void* textureMipSource =
        scene.textureMipLevels.empty()
            ? static_cast<const void*>(&dummyMip)
            : static_cast<const void*>(scene.textureMipLevels.data());
    const void* textureTexelSource =
        scene.textureTexelsRgba8.empty()
            ? static_cast<const void*>(&dummyTexel)
            : static_cast<const void*>(scene.textureTexelsRgba8.data());
    const void* areaLightSource =
        scene.areaLights.empty()
            ? static_cast<const void*>(&dummyAreaLight)
            : static_cast<const void*>(scene.areaLights.data());
    const void* areaLightTriangleSource =
        scene.areaLightTriangles.empty()
            ? static_cast<const void*>(&dummyAreaLightTriangle)
            : static_cast<const void*>(scene.areaLightTriangles.data());
    const void* environmentRowSource =
        scene.environmentRows.empty()
            ? static_cast<const void*>(&dummyEnvironmentRow)
            : static_cast<const void*>(scene.environmentRows.data());
    const void* environmentTexelSource =
        scene.environmentTexels.empty()
            ? static_cast<const void*>(&dummyEnvironmentTexel)
            : static_cast<const void*>(scene.environmentTexels.data());
    const VkDeviceSize packedStagingBytes = std::min(
        std::max(
            {vertexBytes,
             indexBytes,
             triangleMaterialIdBytes,
             materialBytes,
             instanceBytes,
             textureDescriptorBytes,
             textureMipBytes,
             textureTexelBytes,
             areaLightBytes,
             areaLightTriangleBytes,
             environmentRowBytes,
             environmentTexelBytes}),
        kPackedUploadChunkBytes);
    VulkanBuffer packedStaging{
        logicalDevice,
        state_->physical.memoryProperties,
        packedStagingBytes,
        VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
        kHostUploadMemory,
        false};
    const auto uploadPackedBuffer =
        [&](const VulkanBuffer& destination,
            const void* source,
            VkDeviceSize byteCount,
            VkAccessFlags destinationAccess,
            VkPipelineStageFlags destinationStages,
            const char* description) {
            const auto* sourceBytes = static_cast<const std::uint8_t*>(source);
            for (VkDeviceSize offset = 0U; offset < byteCount;) {
                const VkDeviceSize chunkBytes =
                    std::min(byteCount - offset, packedStagingBytes);
                packedStaging.write(
                    sourceBytes + static_cast<std::size_t>(offset), chunkBytes);
                const VkCommandBuffer commands = beginCommands();
                const VkBufferCopy copy{0U, offset, chunkBytes};
                vkCmdCopyBuffer(
                    commands, packedStaging.get(), destination.get(), 1U, &copy);
                const VkBufferMemoryBarrier transferToDestination{
                    VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER,
                    nullptr,
                    VK_ACCESS_TRANSFER_WRITE_BIT,
                    destinationAccess,
                    VK_QUEUE_FAMILY_IGNORED,
                    VK_QUEUE_FAMILY_IGNORED,
                    destination.get(),
                    offset,
                    chunkBytes,
                };
                vkCmdPipelineBarrier(
                    commands,
                    VK_PIPELINE_STAGE_TRANSFER_BIT,
                    destinationStages,
                    0U,
                    0U,
                    nullptr,
                    1U,
                    &transferToDestination,
                    0U,
                    nullptr);
                submitAndWait(description);
                offset += chunkBytes;
            }
        };
    constexpr VkAccessFlags kGeometryDestinationAccess =
        VK_ACCESS_ACCELERATION_STRUCTURE_READ_BIT_KHR | VK_ACCESS_SHADER_READ_BIT;
    constexpr VkPipelineStageFlags kGeometryDestinationStages =
        VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR |
        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT;
    constexpr VkAccessFlags kShaderDestinationAccess = VK_ACCESS_SHADER_READ_BIT;
    constexpr VkPipelineStageFlags kShaderDestinationStage =
        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT;
    uploadPackedBuffer(
        *state_->vertexBuffer,
        scene.vertices.data(),
        vertexBytes,
        kGeometryDestinationAccess,
        kGeometryDestinationStages,
        "packed vertex upload");
    uploadPackedBuffer(
        *state_->indexBuffer,
        scene.triangleIndices.data(),
        indexBytes,
        kGeometryDestinationAccess,
        kGeometryDestinationStages,
        "packed index upload");
    uploadPackedBuffer(
        *state_->triangleMaterialIdBuffer,
        scene.triangleMaterialIds.data(),
        triangleMaterialIdBytes,
        kShaderDestinationAccess,
        kShaderDestinationStage,
        "packed triangle material ID upload");
    uploadPackedBuffer(
        *state_->materialBuffer,
        scene.materials.data(),
        materialBytes,
        kShaderDestinationAccess,
        kShaderDestinationStage,
        "packed material upload");
    uploadPackedBuffer(
        *state_->instanceBuffer,
        &instanceRecord,
        instanceBytes,
        VK_ACCESS_ACCELERATION_STRUCTURE_READ_BIT_KHR,
        VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR,
        "acceleration-structure instance upload");
    uploadPackedBuffer(
        *state_->textureDescriptorBuffer,
        textureDescriptorSource,
        textureDescriptorBytes,
        kShaderDestinationAccess,
        kShaderDestinationStage,
        "packed texture descriptor upload");
    uploadPackedBuffer(
        *state_->textureMipBuffer,
        textureMipSource,
        textureMipBytes,
        kShaderDestinationAccess,
        kShaderDestinationStage,
        "packed texture mip upload");
    uploadPackedBuffer(
        *state_->textureTexelBuffer,
        textureTexelSource,
        textureTexelBytes,
        kShaderDestinationAccess,
        kShaderDestinationStage,
        "packed texture texel upload");
    uploadPackedBuffer(
        *state_->areaLightBuffer,
        areaLightSource,
        areaLightBytes,
        kShaderDestinationAccess,
        kShaderDestinationStage,
        "packed area light upload");
    uploadPackedBuffer(
        *state_->areaLightTriangleBuffer,
        areaLightTriangleSource,
        areaLightTriangleBytes,
        kShaderDestinationAccess,
        kShaderDestinationStage,
        "packed area light triangle upload");
    uploadPackedBuffer(
        *state_->environmentRowBuffer,
        environmentRowSource,
        environmentRowBytes,
        kShaderDestinationAccess,
        kShaderDestinationStage,
        "packed environment row upload");
    uploadPackedBuffer(
        *state_->environmentTexelBuffer,
        environmentTexelSource,
        environmentTexelBytes,
        kShaderDestinationAccess,
        kShaderDestinationStage,
        "packed environment texel upload");
    state_->timings.uploadMilliseconds = elapsedMilliseconds(uploadBegin);

    VkAccelerationStructureGeometryKHR topGeometry{
        VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_GEOMETRY_KHR};
    topGeometry.geometryType = VK_GEOMETRY_TYPE_INSTANCES_KHR;
    topGeometry.geometry.instances.sType =
        VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_GEOMETRY_INSTANCES_DATA_KHR;
    topGeometry.geometry.instances.arrayOfPointers = VK_FALSE;
    topGeometry.geometry.instances.data.deviceAddress =
        state_->instanceBuffer->deviceAddress();
    VkAccelerationStructureBuildGeometryInfoKHR topBuildInfo{
        VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_BUILD_GEOMETRY_INFO_KHR};
    topBuildInfo.type = VK_ACCELERATION_STRUCTURE_TYPE_TOP_LEVEL_KHR;
    topBuildInfo.flags = VK_BUILD_ACCELERATION_STRUCTURE_PREFER_FAST_TRACE_BIT_KHR;
    topBuildInfo.mode = VK_BUILD_ACCELERATION_STRUCTURE_MODE_BUILD_KHR;
    topBuildInfo.geometryCount = 1U;
    topBuildInfo.pGeometries = &topGeometry;
    constexpr std::uint32_t kTopPrimitiveCount = 1U;
    VkAccelerationStructureBuildSizesInfoKHR topSizes{
        VK_STRUCTURE_TYPE_ACCELERATION_STRUCTURE_BUILD_SIZES_INFO_KHR};
    state_->logical.acceleration.getBuildSizes(
        logicalDevice,
        VK_ACCELERATION_STRUCTURE_BUILD_TYPE_DEVICE_KHR,
        &topBuildInfo,
        &kTopPrimitiveCount,
        &topSizes);
    state_->topLevel = std::make_unique<AccelerationStructure>(
        logicalDevice,
        state_->logical.acceleration,
        state_->physical.memoryProperties,
        topSizes.accelerationStructureSize,
        VK_ACCELERATION_STRUCTURE_TYPE_TOP_LEVEL_KHR);
    topBuildInfo.dstAccelerationStructure = state_->topLevel->get();

    const VkDeviceSize scratchRequired =
        std::max(bottomSizes.buildScratchSize, topSizes.buildScratchSize);
    state_->scratchBuffer = std::make_unique<VulkanBuffer>(
        logicalDevice,
        state_->physical.memoryProperties,
        checkedScratchSize(scratchRequired, state_->physical.scratchAlignment),
        VK_BUFFER_USAGE_STORAGE_BUFFER_BIT |
            VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
        true);
    const VkDeviceAddress scratchAddress = alignAddress(
        state_->scratchBuffer->deviceAddress(), state_->physical.scratchAlignment);
    bottomBuildInfo.scratchData.deviceAddress = scratchAddress;
    topBuildInfo.scratchData.deviceAddress = scratchAddress;

    const VulkanClock::time_point buildBegin = VulkanClock::now();
    const VkCommandBuffer buildCommands = beginCommands();
    const VkAccelerationStructureBuildRangeInfoKHR bottomRange{
        primitiveCount,
        0U,
        0U,
        0U,
    };
    const VkAccelerationStructureBuildRangeInfoKHR* bottomRangePointer = &bottomRange;
    state_->logical.acceleration.cmdBuild(
        buildCommands, 1U, &bottomBuildInfo, &bottomRangePointer);
    const VkMemoryBarrier bottomToTopBarrier{
        VK_STRUCTURE_TYPE_MEMORY_BARRIER,
        nullptr,
        VK_ACCESS_ACCELERATION_STRUCTURE_READ_BIT_KHR |
            VK_ACCESS_ACCELERATION_STRUCTURE_WRITE_BIT_KHR,
        VK_ACCESS_ACCELERATION_STRUCTURE_READ_BIT_KHR |
            VK_ACCESS_ACCELERATION_STRUCTURE_WRITE_BIT_KHR,
    };
    vkCmdPipelineBarrier(
        buildCommands,
        VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR,
        VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR,
        0U,
        1U,
        &bottomToTopBarrier,
        0U,
        nullptr,
        0U,
        nullptr);
    const VkAccelerationStructureBuildRangeInfoKHR topRange{
        kTopPrimitiveCount,
        0U,
        0U,
        0U,
    };
    const VkAccelerationStructureBuildRangeInfoKHR* topRangePointer = &topRange;
    state_->logical.acceleration.cmdBuild(
        buildCommands, 1U, &topBuildInfo, &topRangePointer);
    const VkMemoryBarrier buildToTraceBarrier{
        VK_STRUCTURE_TYPE_MEMORY_BARRIER,
        nullptr,
        VK_ACCESS_ACCELERATION_STRUCTURE_WRITE_BIT_KHR,
        VK_ACCESS_ACCELERATION_STRUCTURE_READ_BIT_KHR,
    };
    vkCmdPipelineBarrier(
        buildCommands,
        VK_PIPELINE_STAGE_ACCELERATION_STRUCTURE_BUILD_BIT_KHR,
        VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
        0U,
        1U,
        &buildToTraceBarrier,
        0U,
        nullptr,
        0U,
        nullptr);
    submitAndWait("acceleration-structure build");
    state_->timings.accelerationStructureBuildMilliseconds =
        elapsedMilliseconds(buildBegin);
    // This immutable runtime never updates or rebuilds its acceleration structures.
    state_->scratchBuffer.reset();
    state_->instanceBuffer.reset();
}

VulkanRuntime::VulkanRuntime(
    const PackedSceneData& scene, VulkanRayQueryOptions options) {
    // The backend-neutral contract owns format, reserved-field, finite-value, and ID checks.
    // Device-specific validation below must not drift into a competing scene validator.
    scene.validate();
    implementation_ = std::make_unique<Implementation>(scene, options);
}

VulkanRuntime::~VulkanRuntime() = default;

const std::string& VulkanRuntime::deviceName() const noexcept {
    return implementation_->deviceName();
}

const VulkanRayQuerySetupTimings& VulkanRuntime::setupTimings() const noexcept {
    return implementation_->setupTimings();
}

VulkanValidationReport VulkanRuntime::validationReport() const {
    return implementation_->validationReport();
}

std::size_t VulkanRuntime::triangleCount() const noexcept {
    return implementation_->triangleCount();
}

VkDevice VulkanRuntime::device() const noexcept {
    return implementation_->device();
}

const VkPhysicalDeviceProperties& VulkanRuntime::physicalProperties() const noexcept {
    return implementation_->physicalProperties();
}

const VkPhysicalDeviceMemoryProperties&
VulkanRuntime::memoryProperties() const noexcept {
    return implementation_->memoryProperties();
}

std::uint32_t VulkanRuntime::timestampValidBits() const noexcept {
    return implementation_->timestampValidBits();
}

VkAccelerationStructureKHR
VulkanRuntime::topLevelAccelerationStructure() const noexcept {
    return implementation_->topLevelAccelerationStructure();
}

const VulkanBuffer& VulkanRuntime::vertexBuffer() const noexcept {
    return implementation_->vertexBuffer();
}

const VulkanBuffer& VulkanRuntime::indexBuffer() const noexcept {
    return implementation_->indexBuffer();
}

const VulkanBuffer& VulkanRuntime::triangleMaterialIdBuffer() const noexcept {
    return implementation_->triangleMaterialIdBuffer();
}

const VulkanBuffer& VulkanRuntime::materialBuffer() const noexcept {
    return implementation_->materialBuffer();
}

const VulkanBuffer& VulkanRuntime::textureDescriptorBuffer() const noexcept {
    return implementation_->textureDescriptorBuffer();
}

const VulkanBuffer& VulkanRuntime::textureMipBuffer() const noexcept {
    return implementation_->textureMipBuffer();
}

const VulkanBuffer& VulkanRuntime::textureTexelBuffer() const noexcept {
    return implementation_->textureTexelBuffer();
}

const VulkanBuffer& VulkanRuntime::areaLightBuffer() const noexcept {
    return implementation_->areaLightBuffer();
}

const VulkanBuffer& VulkanRuntime::areaLightTriangleBuffer() const noexcept {
    return implementation_->areaLightTriangleBuffer();
}

const VulkanBuffer& VulkanRuntime::environmentRowBuffer() const noexcept {
    return implementation_->environmentRowBuffer();
}

const VulkanBuffer& VulkanRuntime::environmentTexelBuffer() const noexcept {
    return implementation_->environmentTexelBuffer();
}

std::mutex& VulkanRuntime::operationMutex() noexcept {
    return implementation_->operationMutex();
}

VkCommandBuffer VulkanRuntime::beginCommands() {
    return implementation_->beginCommands();
}

void VulkanRuntime::submitAndWait(const char* description) {
    implementation_->submitAndWait(description);
}

}  // namespace raym0nade::gpu::detail
