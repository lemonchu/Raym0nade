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

struct GbufferData {
    const Face *face;
    vec3 shapeNormal, surfaceNormal, position;
    GbufferData();
};

class Image {
public:
    static const int
            DirectLight = 1,
            IndirectLight = 2,
            DiffuseColor = 4,
            Emission = 8;
    int width, height;
    GbufferData *Gbuffer;
    RadianceData *radiance_d, *radiance_i;
    vec3 *color;
    Image(int width, int height);

    void filter();

    void shade(const Model &model, const vec3 &position, int options);

    void bloom();

    void gammaCorrection();

    void save(const char* file_name);

    ~Image();
};

#endif // IMAGE_H