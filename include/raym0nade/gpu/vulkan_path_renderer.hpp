#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>

#include "raym0nade/gpu/vulkan_ray_query.hpp"
#include "raym0nade/render.hpp"

namespace raym0nade {

class PackedSceneData;

namespace gpu {

struct VulkanPathRenderOptions {
    VulkanRayQueryOptions vulkan;
    std::uint32_t tileWidth{128U};
    std::uint32_t tileHeight{128U};
    std::uint32_t samplesPerBatch{8U};
};

struct VulkanPathRenderTimings {
    // Includes command recording, queue submission/wait, tile readback, and Film assembly.
    double hostRenderMilliseconds{0.0};
    // Sum of per-queue device timestamps enclosing compute work for every tile. With multiple
    // queues this is aggregate queue-busy time, not elapsed GPU wall-clock time.
    double gpuDispatchMilliseconds{0.0};
    bool gpuTimestampAvailable{false};
    std::uint64_t dispatchCount{0U};
    std::uint32_t computeQueueCount{1U};
};

struct VulkanPathRenderResult {
    Film film;
    RenderStats stats;
    VulkanPathRenderTimings timings;
};

// Persistent, headless Vulkan path renderer for one immutable packed scene. A render is
// subdivided into bounded tiles and deterministic sample batches. The returned Film is ready
// for the same export and post-processing path used by the CPU reference renderer.
class VulkanPathRenderer {
public:
    VulkanPathRenderer(
        const PackedSceneData& scene,
        const std::filesystem::path& spirvPath,
        VulkanPathRenderOptions options = {});
    ~VulkanPathRenderer();

    VulkanPathRenderer(const VulkanPathRenderer&) = delete;
    VulkanPathRenderer& operator=(const VulkanPathRenderer&) = delete;
    VulkanPathRenderer(VulkanPathRenderer&&) = delete;
    VulkanPathRenderer& operator=(VulkanPathRenderer&&) = delete;

    [[nodiscard]] const std::string& deviceName() const noexcept;
    [[nodiscard]] std::uint32_t computeQueueCount() const noexcept;
    [[nodiscard]] const VulkanRayQuerySetupTimings& setupTimings() const noexcept;
    [[nodiscard]] VulkanValidationReport validationReport() const;
    [[nodiscard]] VulkanPathRenderResult render(const RenderSettings& settings);

private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};

}  // namespace gpu
}  // namespace raym0nade
