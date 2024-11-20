#ifndef IMAGE_H
#define IMAGE_H

#include "geometry.h"

struct PixelData {
    vec3 color;
    float depth;
    PixelData();
    PixelData(const vec3 &color, const float &depth);
};

PixelData operator + (const PixelData &A, const PixelData &B);

class Image {
public:
    unsigned int width, height;
    PixelData* buffer;

    Image(unsigned int width, unsigned int height);

    void save(const char* file_name);

    void normalize();

    ~Image();
};

#endif // IMAGE_H