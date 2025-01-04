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

class Photo {
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
        DoDepthFieldBlur = 1024
    };
    int width, height;
    float exposure, focus, CoC;
    vec3 cameraPosition;
    HitInfo *Gbuffer;
    RadianceData *radiance_Dd, *radiance_Ds, *radiance_Id, *radiance_Is;
    vec3 *pixelarray;
    Photo(int width, int height);
    Photo(const char* file_name);

    void filter();

    void shade(int options);

    void bloom();

    void depthFeildBlur();

    void FXAA();

    void gammaCorrection();
    void reverseGammaCorrection();

    void save(const char* file_name);

    void postProcessing(int shadeOptions);

    ~Photo();

private:
    void load(const char* file_name);
};

void accumulateInwardRadiance(const vec3 &baseColor, const LightSample &sample,
                              RadianceData &radiance_d, RadianceData &radiance_s);

#endif // IMAGE_H