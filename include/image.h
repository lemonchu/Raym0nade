#ifndef IMAGE_H
#define IMAGE_H

#include "geometry.h"
#include "component.h"

struct RadianceData {
    vec3 radiance;
    float Var;
    int sampleCount;
    RadianceData();
};

struct GbufferData {
    vec3 diffuseColor, shapeNormal, surfaceNormal, position, emission;
    const Face *face;
    float depth;
    GbufferData();
};

class Image {
public:
    unsigned int width, height;
    GbufferData *Gbuffer;
    RadianceData *radiance_d, *radiance_i;
    vec3 *color;
    Image(unsigned int width, unsigned int height);

    void normalizeRadiance();

    void filter();

    void shade();

    void bloom();

    void gammaCorrection();

    void save(const char* file_name);

    ~Image();
};

#endif // IMAGE_H