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
};

class Image {
public:
    int width, height;
    PixelData* buffer;

    Image(int width, int height) : width(width), height(height) {
        buffer = new PixelData[width * height];
    }

    int getWidth() const { return width; }
    int getHeight() const { return height; }

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
};

#endif // IMAGE_H