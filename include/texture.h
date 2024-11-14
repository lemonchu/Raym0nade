#ifndef TEXTURE_H
#define TEXTURE_H

#include <cstdint>
#include <cstdlib>
#include <png.h>
#include <stdexcept>
#include <cstring>
#include <vector>
#include <iostream>

#define MAX_IMAGES 30

class Texture {
public:
    int width;
    int height;
    std::vector<uint8_t> images[MAX_IMAGES];
    bool enabled[MAX_IMAGES]; // Array to indicate if each image is enabled

    Texture() : width(0), height(0) {
        memset(enabled, 0, sizeof(enabled));
    }

    int isEnabled(int index) {
        return enabled[index];
    }

    void addImage(int index, const std::vector<uint8_t>& imageData) {
        if (index < 0 || index >= MAX_IMAGES) {
            throw std::out_of_range("Index out of range");
        }
        if (imageData.size() != width * height * 3) {
            throw std::invalid_argument("Image data size does not match texture dimensions");
        }
        images[index] = imageData;
        enabled[index] = true; // Enable the image when added
    }

    [[nodiscard]] const std::vector<uint8_t>& getImage(int index) const {
        if (index < 0 || index >= MAX_IMAGES) {
            throw std::out_of_range("Index out of range");
        }
        if (!enabled[index]) {
            throw std::runtime_error("Image is not enabled");
        }
        return images[index];
    }

    void loadImageFormFile(int index, const std::string& filename) {
        if (index < 0 || index >= MAX_IMAGES) {
            throw std::out_of_range("Index out of range");
        }

        FILE *fp = fopen(filename.c_str(), "rb");
        if (!fp) throw std::runtime_error("Failed to open file for reading");

        png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        if (!png) throw std::runtime_error("Failed to create PNG read struct");

        png_infop info = png_create_info_struct(png);
        if (!info) throw std::runtime_error("Failed to create PNG info struct");

        if (setjmp(png_jmpbuf(png))) throw std::runtime_error("Error during PNG reading");

        png_init_io(png, fp);
        png_read_info(png, info);

        int width = png_get_image_width(png, info);
        int height = png_get_image_height(png, info);
        png_byte color_type = png_get_color_type(png, info);
        png_byte bit_depth = png_get_bit_depth(png, info);

        if (this->width == 0) {
            // If this is the first image, set the dimensions
            this->width = width;
            this->height = height;
        }

        if (width != this->width || height != this->height) {
            throw std::runtime_error("Image dimensions do not match texture requirements");
        }

        if (color_type == PNG_COLOR_TYPE_GRAY) {
            png_set_gray_to_rgb(png);
        } else if (color_type != PNG_COLOR_TYPE_RGB) {
            throw std::runtime_error("Unsupported image format");
        }

        png_read_update_info(png, info);

        std::vector<uint8_t> imageData(width * height * 3);
        for (int y = 0; y < height; y++) {
            png_read_row(png, &imageData[y * width * 3], nullptr);
        }

        fclose(fp);
        png_destroy_read_struct(&png, &info, nullptr);

        images[index].resize(width * height * 3);
        images[index] = imageData;
        enabled[index] = true; // Enable the image when loaded
    }

    void saveImageToFile(int index, const std::string& filename) const {
        if (index < 0 || index >= MAX_IMAGES) {
            throw std::out_of_range("Index out of range");
        }
        if (!enabled[index]) {
            throw std::runtime_error("Image is not enabled");
        }

        FILE *fp = fopen(filename.c_str(), "wb");
        if (!fp) throw std::runtime_error("Failed to open file for writing");

        png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        if (!png) throw std::runtime_error("Failed to create PNG write struct");

        png_infop info = png_create_info_struct(png);
        if (!info) throw std::runtime_error("Failed to create PNG info struct");

        if (setjmp(png_jmpbuf(png))) throw std::runtime_error("Error during PNG creation");

        png_init_io(png, fp);

        png_set_IHDR(
                png,
                info,
                width, height,
                8, PNG_COLOR_TYPE_RGB,
                PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_DEFAULT,
                PNG_FILTER_TYPE_DEFAULT
        );
        png_write_info(png, info);

        for (int y = 0; y < height; y++)
            png_write_row(png, &images[index][y * width * 3]);

        png_write_end(png, nullptr);

        fclose(fp);
        png_destroy_write_struct(&png, &info);
    }
};

#endif // TEXTURE_H