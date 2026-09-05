#pragma once

#include <vulkan/vulkan.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "raym0nade/gpu/vulkan_ray_query.hpp"

namespace raym0nade {

class PackedSceneData;

namespace gpu::detail {

using VulkanClock = std::chrono::steady_clock;

[[nodiscard]] double elapsedMilliseconds(VulkanClock::time_point begin) noexcept;
void requireVulkanSuccess(VkResult result, const char* operation);
[[nodiscard]] std::vector<std::uint32_t> readSpirvFile(
    const std::filesystem::path& path);
[[nodiscard]] VkDeviceSize checkedVulkanByteSize(
    std::size_t count, std::size_t stride, const char* description);

template <typename Handle>
class UniqueVulkanHandle {
public:
    UniqueVulkanHandle() = default;

    UniqueVulkanHandle(Handle handle, std::function<void(Handle)> deleter)
        : handle_(handle), deleter_(std::move(deleter)) {}

    UniqueVulkanHandle(const UniqueVulkanHandle&) = delete;
    UniqueVulkanHandle& operator=(const UniqueVulkanHandle&) = delete;

    UniqueVulkanHandle(UniqueVulkanHandle&& other) noexcept
        : handle_(std::exchange(other.handle_, Handle{})),
          deleter_(std::move(other.deleter_)) {}

    UniqueVulkanHandle& operator=(UniqueVulkanHandle&& other) noexcept {
        if (this != &other) {
            reset();
            handle_ = std::exchange(other.handle_, Handle{});
            deleter_ = std::move(other.deleter_);
        }
        return *this;
    }

    ~UniqueVulkanHandle() {
        reset();
    }

    void reset(Handle handle = Handle{}, std::function<void(Handle)> deleter = {}) {
        if (handle_ != Handle{} && deleter_) {
            deleter_(handle_);
        }
        handle_ = handle;
        deleter_ = std::move(deleter);
    }

    [[nodiscard]] Handle get() const noexcept {
        return handle_;
    }

private:
    Handle handle_{};
    std::function<void(Handle)> deleter_;
};

class VulkanBuffer {
public:
    VulkanBuffer(
        VkDevice device,
        const VkPhysicalDeviceMemoryProperties& memoryProperties,
        VkDeviceSize size,
        VkBufferUsageFlags usage,
        VkMemoryPropertyFlags requiredMemoryProperties,
        bool needsDeviceAddress);
    ~VulkanBuffer();

    VulkanBuffer(const VulkanBuffer&) = delete;
    VulkanBuffer& operator=(const VulkanBuffer&) = delete;
    VulkanBuffer(VulkanBuffer&&) = delete;
    VulkanBuffer& operator=(VulkanBuffer&&) = delete;

    void write(const void* source, VkDeviceSize byteCount) const;
    void read(void* destination, VkDeviceSize byteCount) const;

    [[nodiscard]] VkBuffer get() const noexcept;
    [[nodiscard]] VkDeviceSize size() const noexcept;
    [[nodiscard]] VkDeviceAddress deviceAddress() const;

private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};

// Owns one Vulkan device, one immutable packed-geometry upload, and its persistent BLAS/TLAS.
// Operations share one command buffer and queue, so callers must hold operationMutex() from
// beginCommands() through submitAndWait(). This private boundary is reusable by later renderers
// without exposing Vulkan types through the public include tree.
class VulkanRuntime {
public:
    VulkanRuntime(const PackedSceneData& scene, VulkanRayQueryOptions options);
    ~VulkanRuntime();

    VulkanRuntime(const VulkanRuntime&) = delete;
    VulkanRuntime& operator=(const VulkanRuntime&) = delete;
    VulkanRuntime(VulkanRuntime&&) = delete;
    VulkanRuntime& operator=(VulkanRuntime&&) = delete;

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

    [[nodiscard]] std::mutex& operationMutex() noexcept;
    [[nodiscard]] VkCommandBuffer beginCommands();
    void submitAndWait(const char* description);

private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};

}  // namespace gpu::detail
}  // namespace raym0nade
