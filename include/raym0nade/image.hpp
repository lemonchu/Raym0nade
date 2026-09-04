#pragma once

#include <filesystem>
#include <vector>

#include "raym0nade/model.hpp"
#include "raym0nade/sampling.hpp"

namespace raym0nade {

struct RadianceData {
    vec3 radiance{0.0F};
    float varianceAccumulator{0.0F};
};

class Film {
public:
    enum ShadeOption {
        baseColor = 1,
        emission = 2,
        directLight = 4,
        indirectLight = 8,
        diffuse = 16,
        specular = 32,
        shapeNormal = 64,
        surfaceNormal = 128,
        directDiffuse = directLight | diffuse,
        directSpecular = directLight | specular,
        indirectDiffuse = indirectLight | diffuse,
        indirectSpecular = indirectLight | specular,
        full = directLight | indirectLight | diffuse | specular | baseColor | emission,
        bloomEnabled = 256,
        fxaaEnabled = 512,
        depthOfFieldEnabled = 1024,
    };

    Film(int width, int height);
    explicit Film(const std::filesystem::path& filename);

    Film(const Film&) = delete;
    Film& operator=(const Film&) = delete;
    Film(Film&&) noexcept = default;
    Film& operator=(Film&&) noexcept = default;
    ~Film() = default;

    void spatialClamp();
    void filter();
    void shade(int options);
    void bloom();
    void depthOfFieldBlur();
    void applyFxaa();
    void gammaCorrect();
    void reverseGammaCorrect();
    void postProcess(int shadeOptions);
    void save(const std::filesystem::path& filename) const;

    [[nodiscard]] int width() const noexcept;
    [[nodiscard]] int height() const noexcept;

    float exposure{1.0F};
    float focusDistance{0.0F};
    float circleOfConfusion{0.0F};
    vec3 cameraPosition{0.0F};
    std::vector<HitInfo> gBuffer;
    std::vector<RadianceData> directDiffuseRadiance;
    std::vector<RadianceData> directSpecularRadiance;
    std::vector<RadianceData> indirectDiffuseRadiance;
    std::vector<RadianceData> indirectSpecularRadiance;
    std::vector<vec3> pixels;

private:
    void load(const std::filesystem::path& filename);

    int width_{0};
    int height_{0};
};

void accumulateInwardRadiance(
    const vec3& baseColor,
    const LightSample& sample,
    RadianceData& diffuseRadiance,
    RadianceData& specularRadiance);

}  // namespace raym0nade
