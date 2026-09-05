#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include "raym0nade/gpu/vulkan_ray_query.hpp"
#include "raym0nade/render_contract.hpp"

namespace raym0nade {

class PackedSceneData;

namespace gpu {

struct VulkanPrimaryRenderTimings {
    // Includes queue submission/wait and host readback. This is a host wall-clock measurement.
    double dispatchAndReadbackMilliseconds{0.0};
    // Set only when the selected queue supports valid Vulkan timestamp queries.
    double gpuDispatchMilliseconds{0.0};
    bool gpuTimestampAvailable{false};
};

struct VulkanPrimaryRenderResult {
    LinearImage image;
    VulkanPrimaryRenderTimings timings;
};

// Persistent, headless AMD Vulkan primary-AOV renderer for one immutable packed scene.
// Construction uploads scene buffers and builds one BLAS plus one TLAS. render() creates
// primary rays on the device and returns a row-major linear RGB image.
// ShapeNormal supports every scene accepted by this backend. BaseColor currently requires
// every referenced material to be opaque and to have no diffuse texture; render() throws
// std::invalid_argument when that restricted BaseColor contract is not satisfied.
class VulkanPrimaryRenderer {
public:
    VulkanPrimaryRenderer(
        const PackedSceneData& scene,
        const std::filesystem::path& spirvPath,
        VulkanRayQueryOptions options = {});
    ~VulkanPrimaryRenderer();

    VulkanPrimaryRenderer(const VulkanPrimaryRenderer&) = delete;
    VulkanPrimaryRenderer& operator=(const VulkanPrimaryRenderer&) = delete;
    VulkanPrimaryRenderer(VulkanPrimaryRenderer&&) = delete;
    VulkanPrimaryRenderer& operator=(VulkanPrimaryRenderer&&) = delete;

    [[nodiscard]] const std::string& deviceName() const noexcept;
    [[nodiscard]] const VulkanRayQuerySetupTimings& setupTimings() const noexcept;
    [[nodiscard]] VulkanValidationReport validationReport() const;
    [[nodiscard]] VulkanPrimaryRenderResult render(const PrimaryRenderRequest& request);

private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};

}  // namespace gpu
}  // namespace raym0nade
