#pragma once

#include <cstdint>
#include <filesystem>

#include "raym0nade/model.hpp"
#include "raym0nade/render_contract.hpp"

namespace raym0nade {

struct RenderSettings {
    vec3 position{0.0F};
    vec3 direction{0.0F, 0.0F, 1.0F};
    vec3 up{0.0F, 1.0F, 0.0F};
    vec3 right{1.0F, 0.0F, 0.0F};
    float pixelScale{1.0e-3F};
    float focusDistance{0.0F};
    float circleOfConfusion{0.0F};
    float exposure{1.0F};
    float directLightProbability{0.7F};
    int width{640};
    int height{360};
    int samplesPerPixel{1};
    int threadCount{0};
    std::uint32_t seed{0};
    std::filesystem::path outputPrefix{"output/render"};

    void validate() const;
    [[nodiscard]] int resolvedThreadCount() const noexcept;
};

struct RenderStats {
    double renderSeconds{0.0};
    double totalSeconds{0.0};
    std::uint64_t directLightSamples{0};
};

struct CpuPrimaryRenderOptions {
    // Zero selects hardware concurrency; positive values request an explicit worker count.
    int threadCount{1};

    void validate() const;
    [[nodiscard]] std::uint32_t resolvedThreadCount(std::uint32_t imageHeight) const noexcept;
};

[[nodiscard]] RenderStats renderToFiles(const Model& model, const RenderSettings& settings);
[[nodiscard]] LinearImage renderPrimaryAovCpu(
    const Model& model, const PrimaryRenderRequest& request);
[[nodiscard]] LinearImage renderPrimaryAovCpu(
    const Model& model,
    const PrimaryRenderRequest& request,
    const CpuPrimaryRenderOptions& options);

}  // namespace raym0nade
