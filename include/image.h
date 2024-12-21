#ifndef IMAGE_H
#define IMAGE_H

#include "geometry.h"
#include "component.h"
#include "model.h"

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
        Full = DirectLight | IndirectLight | Diffuse | Specular | BaseColor | Emission
    };
    int width, height;
    HitInfo *Gbuffer;
    RadianceData *radiance_Dd, *radiance_Ds, *radiance_Id, *radiance_Is;
    vec3 *color;
    Image(int width, int height);

    void filter();

    void shade(float exposure, int options);

    void FXAA();

    void bloom();

    void gammaCorrection();

    void save(const char* file_name);

    ~Image();
};

#endif // IMAGE_H