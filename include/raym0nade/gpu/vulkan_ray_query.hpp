#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace raym0nade {

class PackedSceneData;

namespace gpu {

inline constexpr std::uint32_t kInvalidRayQueryPrimitiveId = 0xffffffffU;

// Fixed host/shader ABI. The fourth component stores the valid parametric interval.
struct alignas(16) VulkanRayQueryRay {
    std::array<float, 4> originAndTMin{};
    std::array<float, 4> directionAndTMax{};
};

// Fixed host/shader ABI. `hit` is an integer so the record has no C++ bool ABI dependency.
struct alignas(16) VulkanRayQueryHit {
    std::uint32_t hit{0U};
    std::uint32_t primitiveId{kInvalidRayQueryPrimitiveId};
    float distance{0.0F};
    float barycentricU{0.0F};
    float barycentricV{0.0F};
    std::array<std::uint32_t, 3> reserved{};

    [[nodiscard]] bool hasHit() const noexcept {
        return hit != 0U;
    }
};

struct VulkanRayQuerySetupTimings {
    // Host wall-clock measurements. These are bring-up diagnostics, not GPU timestamps.
    double uploadMilliseconds{0.0};
    double accelerationStructureBuildMilliseconds{0.0};
};

struct VulkanRayQueryBatch {
    std::vector<VulkanRayQueryHit> hits;
    // Includes host upload, queue submission/wait, dispatch, and host readback.
    // This is a host wall-clock measurement, not a GPU timestamp.
    double dispatchAndReadbackMilliseconds{0.0};
};

struct VulkanValidationReport {
    bool requested{false};
    bool enabled{false};
    bool synchronizationValidationEnabled{false};
    std::uint32_t errorCount{0U};
    std::uint32_t warningCount{0U};
    std::vector<std::string> messages;
};

struct VulkanRayQueryOptions {
    bool requestValidation{false};
    std::uint64_t fenceTimeoutNanoseconds{60'000'000'000ULL};
    // Upper bound for one packed RGBA8 texel page. The runtime rounds down to a power of two
    // and also respects maxStorageBufferRange. Small values are useful for paging tests.
    std::uint64_t textureTexelPageBytes{256ULL * 1024ULL * 1024ULL};
};

// Persistent, headless AMD Vulkan Ray Query intersector for one immutable packed scene.
// Construction uploads geometry and builds one BLAS plus one TLAS. Each intersect() call
// submits the complete batch once; it never synchronizes per ray. This low-level intersector
// rejects alpha-cutout scenes because its compact result ABI has no material sampling contract.
class VulkanRayQueryIntersector {
public:
    VulkanRayQueryIntersector(
        const PackedSceneData& scene,
        const std::filesystem::path& spirvPath,
        VulkanRayQueryOptions options = {});
    ~VulkanRayQueryIntersector();

    VulkanRayQueryIntersector(const VulkanRayQueryIntersector&) = delete;
    VulkanRayQueryIntersector& operator=(const VulkanRayQueryIntersector&) = delete;
    VulkanRayQueryIntersector(VulkanRayQueryIntersector&&) = delete;
    VulkanRayQueryIntersector& operator=(VulkanRayQueryIntersector&&) = delete;

    [[nodiscard]] const std::string& deviceName() const noexcept;
    [[nodiscard]] const VulkanRayQuerySetupTimings& setupTimings() const noexcept;
    [[nodiscard]] VulkanValidationReport validationReport() const;
    [[nodiscard]] VulkanRayQueryBatch intersect(
        const std::vector<VulkanRayQueryRay>& rays);

private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};

static_assert(sizeof(VulkanRayQueryRay) == 32U, "Vulkan ray input must match the shader ABI.");
static_assert(alignof(VulkanRayQueryRay) == 16U);
static_assert(offsetof(VulkanRayQueryRay, originAndTMin) == 0U);
static_assert(offsetof(VulkanRayQueryRay, directionAndTMax) == 16U);
static_assert(std::is_standard_layout_v<VulkanRayQueryRay>);
static_assert(std::is_trivially_copyable_v<VulkanRayQueryRay>);

static_assert(sizeof(VulkanRayQueryHit) == 32U, "Vulkan ray result must match the shader ABI.");
static_assert(alignof(VulkanRayQueryHit) == 16U);
static_assert(offsetof(VulkanRayQueryHit, hit) == 0U);
static_assert(offsetof(VulkanRayQueryHit, primitiveId) == 4U);
static_assert(offsetof(VulkanRayQueryHit, distance) == 8U);
static_assert(offsetof(VulkanRayQueryHit, barycentricU) == 12U);
static_assert(offsetof(VulkanRayQueryHit, barycentricV) == 16U);
static_assert(offsetof(VulkanRayQueryHit, reserved) == 20U);
static_assert(std::is_standard_layout_v<VulkanRayQueryHit>);
static_assert(std::is_trivially_copyable_v<VulkanRayQueryHit>);

}  // namespace gpu
}  // namespace raym0nade
