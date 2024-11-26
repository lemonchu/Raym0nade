#include <vector>
#include <png.h>
#include <iostream>
#include "image.h"

PixelData::PixelData() {
    color = vec4(0);
    depth = 0;
}

Image::Image(unsigned int width, unsigned int height) : buffer(nullptr), width(width), height(height) {
    buffer = new PixelData[width * height];
}

void Image::save(const char *file_name) {

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);

    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        std::cerr << "Error during png creation" << std::endl;
        return;
    }

    FILE *fp;
    if ((fp = fopen(file_name, "wb")) == nullptr) {
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
    for (unsigned int y = 0; y < height; ++y)
        for (unsigned int x = 0; x < width; ++x) {
            int id = y * width + x;
            png_byte *row = &image_data[id * 3];
            buffer[id].color /= buffer[id].color[3];
            for (int k = 0; k < 3; ++k) {
                float buf = std::min(std::max(buffer[id].color[k], 0.0f), 1.0f);
                buf = std::pow(buf,1.0f / gamma);
                row[k] = static_cast<png_byte>(buf * 255);
            }

        }

    for (unsigned int y = 0; y < height; ++y)
        png_write_row(png_ptr, &image_data[y * width * 3]);

    png_write_end(png_ptr, nullptr);

    fclose(fp);
    png_destroy_write_struct(&png_ptr, &info_ptr);
}

Image::~Image() {
    delete[] buffer;
}