#ifndef IMAGE_H
#define IMAGE_H

#include "geometry.h"

struct PixelData {
    vec4 color;
    float depth;
    PixelData();
};

PixelData operator + (const PixelData &A, const PixelData &B);

class Image {
public:
    unsigned int width, height;
    PixelData* buffer;

    Image(unsigned int width, unsigned int height);

    void save(const char* file_name);

    ~Image();
};

#endif // IMAGE_H