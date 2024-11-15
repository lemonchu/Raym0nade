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

        std::cout << "Loading image from file: " << filename << std::endl;

        FILE *fp = fopen(filename.c_str(), "rb");
        if (!fp) throw std::runtime_error("Failed to open file for reading");

        png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        if (!png) throw std::runtime_error("Failed to create PNG read struct");

        png_infop info = png_create_info_struct(png);
        if (!info) throw std::runtime_error("Failed to create PNG info struct");

        if (setjmp(png_jmpbuf(png))) throw std::runtime_error("Error during PNG reading");

        png_init_io(png, fp);
        png_read_info(png, info);

        int imgWidth = png_get_image_width(png, info);
        int imgHeight = png_get_image_height(png, info);
        png_byte color_type = png_get_color_type(png, info);
        png_byte bit_depth = png_get_bit_depth(png, info);

        if (bit_depth == 16) {
            png_set_strip_16(png);
        }

        if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) {
            png_set_expand_gray_1_2_4_to_8(png);
        }

        if (color_type == PNG_COLOR_TYPE_PALETTE) {
            png_set_palette_to_rgb(png);
        }

        if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth >= 8) {
            png_set_gray_to_rgb(png);
        }

        if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
            png_set_gray_to_rgb(png);
        }

        if (color_type == PNG_COLOR_TYPE_RGBA) {
            png_set_strip_alpha(png);
        }

        png_read_update_info(png, info);

        std::vector<png_byte> row(imgWidth * 3); // 3 bytes per pixel (RGB)
        std::vector<uint8_t> imageData(imgWidth * imgHeight * 3);

        for (int y = 0; y < imgHeight; y++) {
            png_read_row(png, row.data(), nullptr);
            for (int x = 0; x < imgWidth * 3; x++) {
                imageData[y * imgWidth * 3 + x] = row[x];
            }
        }

        fclose(fp);
        png_destroy_read_struct(&png, &info, nullptr);

        if (width == 0 && height == 0) {
            width = imgWidth;
            height = imgHeight;
        } else if (width != imgWidth || height != imgHeight) {
            throw std::runtime_error("Image dimensions do not match texture requirements");
        }

        images[index] = imageData;
        enabled[index] = true; // Enable the image when loaded
    }
};

#endif // TEXTURE_H