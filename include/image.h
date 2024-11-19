#ifndef IMAGE_H
#define IMAGE_H

#include <glm/glm.hpp>
#include <png.h>
#include <iostream>

struct PixelData {
    glm::vec<3, float> color;
    float depth;
    PixelData() {
        color = glm::vec<3, float>(0);
        depth = 0;
    }
    PixelData(const glm::vec<3, float> &color, float depth) : color(color), depth(depth) {}
};

PixelData operator + (const PixelData &A, const PixelData &B) {
    return PixelData(
        A.color + B.color,
        A.depth + B.depth
    );
}

class Image {
public:
    unsigned int width, height;
    PixelData* buffer;

    Image(unsigned int width, unsigned int height) : buffer(nullptr), width(width), height(height) {
        buffer = new PixelData[width * height];
    }

    void save(const char* file_name) {

        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        png_infop info_ptr = png_create_info_struct(png_ptr);

        if (setjmp(png_jmpbuf(png_ptr))) {
            png_destroy_write_struct(&png_ptr, &info_ptr);
            std::cerr << "Error during png creation" << std::endl;
            return;
        }

        FILE* fp = fopen(file_name, "wb");
        if (!fp) {
            png_destroy_write_struct(&png_ptr, &info_ptr);
            std::cerr << "Could not open file for writing" << std::endl;
            return;
        }

        png_init_io(png_ptr, fp);

        png_set_IHDR(
                png_ptr,
                info_ptr,
                width,
                height,
                8,
                PNG_COLOR_TYPE_RGB,
                PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_DEFAULT,
                PNG_FILTER_TYPE_DEFAULT
        );

        png_write_info(png_ptr, info_ptr);

        std::vector<png_byte> image_data;
        image_data.resize(width * height * 3);
        for (int y = 0; y < height; ++y)
            for (int x = 0; x < width; ++x){
                int id = y * width + x;
                png_byte* row = &image_data[id * 3];
                row[0] = buffer[id].color[0] * 255;
                row[1] = buffer[id].color[1] * 255;
                row[2] = buffer[id].color[2] * 255;
            }

        for (int y = 0; y < height; ++y)
            png_write_row(png_ptr, &image_data[y * width * 3]);


        png_write_end(png_ptr, nullptr);

        fclose(fp);
        png_destroy_write_struct(&png_ptr, &info_ptr);
    }

    void clear() {
        for (int i = 0; i < width * height; ++i)
            buffer[i] = PixelData();
    }

    void normalize() {
        float maxColor = -1e9 , minColor = 1e9;
        for (int i = 0; i < width * height; ++i) {
            if (buffer[i].depth < 1e-6) continue;
            maxColor = std::max(maxColor, std::max(buffer[i].color[0], std::max(buffer[i].color[1], buffer[i].color[2])));
            minColor = std::min(minColor, std::min(buffer[i].color[0], std::min(buffer[i].color[1], buffer[i].color[2])));
        }
        for (int i = 0; i < width * height; ++i) {
            if (buffer[i].depth < 1e-6) continue;
            buffer[i].color = (buffer[i].color - glm::vec<3, float>(minColor)) / (maxColor - minColor);
        }

    }
};

#endif // IMAGE_H