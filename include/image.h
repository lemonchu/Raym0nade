#ifndef IMAGE_H
#define IMAGE_H

#include "geometry.h"
#include "component.h"
#include "model.h"
#include "sampling.h"

struct RadianceData {
    vec3 radiance;
    float Var;
    RadianceData();
};

class Image {
public:
    enum ShadeOption {
        BaseColor = 1,
        Emission = 2,
        DirectLight = 4,
        IndirectLight = 8,
        Diffuse = 16,
        Specular = 32,
        shapeNormal = 64,
        surfaceNormal = 128,
        Direct_Diffuse = DirectLight | Diffuse,
        Direct_Specular = DirectLight | Specular,
        Indirect_Diffuse = IndirectLight | Diffuse,
        Indirect_Specular = IndirectLight | Specular,
        Full = DirectLight | IndirectLight | Diffuse | Specular | BaseColor | Emission,
        DoBloom = 256,
        DoFXAA = 512,
    };
    int width, height;
    HitInfo *Gbuffer;
    RadianceData *radiance_Dd, *radiance_Ds, *radiance_Id, *radiance_Is;
    vec3 *pixelarray;
    Image(int width, int height);
    Image(const char* file_name);

    void filter();

    void shade(float exposure, int options);

    void FXAA();

    void bloom();

    void gammaCorrection();
    void reverseGammaCorrection();

    void save(const char* file_name);

    void postProcessing(int shadeOptions, float exposure);

    ~Image();

private:
    void load(const char* file_name);
};

void accumulateInwardRadiance(const vec3 &baseColor, const LightSample &sample,
                              RadianceData &radiance_d, RadianceData &radiance_s);

#endif // IMAGE_H