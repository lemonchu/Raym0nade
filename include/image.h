#ifndef IMAGE_H
#define IMAGE_H

#include "geometry.h"
#include "component.h"

struct PixelData {
    vec3 inradiance, diffuseColor, shapeNormal, surfaceNormal, position, emission;
    const Face *face;
    int sampleCount;
    float depth, Epower2, Var;
    PixelData();
};

PixelData operator + (const PixelData &A, const PixelData &B);

class Image {
private:
    void filterVar();
    void filterInradiance(int step);
public:
    unsigned int width, height;
    PixelData* buffer;

    Image(unsigned int width, unsigned int height);

    void filterInradiance();

    void save(const char* file_name);

    ~Image();
};

const float gamma = 2.2f;

#endif // IMAGE_H