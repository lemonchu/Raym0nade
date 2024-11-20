#include <assimp/material.h>
#include <cstdint>
#include <cstdlib>
#include <png.h>
#include <turbojpeg.h>
#include <stdexcept>
#include <cstring>
#include <vector>
#include <iostream>
#include "texture.h"

Material::Material() : width(0), height(0) {
    memset(enabled, 0, sizeof(enabled));
}

[[nodiscard]] const std::vector<uint8_t>& Material::getImage(int index) const {
    if (index < 0 || index >= AI_TEXTURE_TYPE_MAX) {
        throw std::out_of_range("Index out of range");
    }
    if (!enabled[index]) {
        throw std::runtime_error("Image is not enabled");
    }
    return texture[index];
}

bool Material::loadImageFromPNG(std::vector<uint8_t> &imageData, const std::string& filename) {

    std::cout << "Loading image from file: " << filename << std::endl;

    FILE *fp;
    if ((fp = fopen(filename.c_str(), "rb")) == nullptr)
        throw std::runtime_error("Failed to open file for reading");

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
    imageData.resize(imgWidth * imgHeight * 3);

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

    return true;
}

bool Material::loadImageFromJPG(std::vector<uint8_t> &imageData, const std::string& filename) {
    // Open JPEG file
    FILE* jpegFile = fopen(filename.c_str(), "rb");
    if (jpegFile == nullptr) {
        std::cerr << "Error opening JPEG file: " << filename << std::endl;
        return false;
    }

    // Create TurboJPEG decompressor instance
    tjhandle tjInstance = tjInitDecompress();
    if (tjInstance == nullptr) {
        std::cerr << "Error initializing TurboJPEG decompressor" << std::endl;
        fclose(jpegFile);
        return false;
    }

    // Read JPEG file content into buffer
    fseek(jpegFile, 0, SEEK_END);
    unsigned long jpegSize = ftell(jpegFile);
    fseek(jpegFile, 0, SEEK_SET);
    unsigned char* jpegBuf = (unsigned char*)malloc(jpegSize);
    if (jpegBuf == nullptr) {
        std::cerr << "Memory allocation failure" << std::endl;
        tjDestroy(tjInstance);
        fclose(jpegFile);
        return false;
    }
    fread(jpegBuf, 1, jpegSize, jpegFile);
    fclose(jpegFile); // Close file as content is already read into buffer

    // Get JPEG image dimensions and components
    int imgWidth, imgHeight, jpegSubsamp;
    if (tjDecompressHeader2(tjInstance, jpegBuf, jpegSize, &imgWidth, &imgHeight, &jpegSubsamp) < 0) {
        std::cerr << "Error getting JPEG image header" << std::endl;
        free(jpegBuf);
        tjDestroy(tjInstance);
        return false;
    }

    // Set output buffer size based on image dimensions and components
    int pixelSize = 3; // RGB has 3 components
    int pitch = imgWidth * pixelSize;
    imageData.resize(imgHeight * pitch);

    // Decompress JPEG image to RGB format
    if (tjDecompress2(tjInstance, jpegBuf, jpegSize, &imageData[0], imgWidth, pitch, imgHeight, TJPF_RGB, TJFLAG_FASTDCT) < 0) {
        std::cerr << "Error decompressing JPEG image" << std::endl;
        free(jpegBuf);
        tjDestroy(tjInstance);
        return false;
    }

    // Free resources
    free(jpegBuf);
    tjDestroy(tjInstance);

    if (width == 0 && height == 0) {
        width = imgWidth;
        height = imgHeight;
    } else if (width != imgWidth || height != imgHeight) {
        throw std::runtime_error("Image dimensions do not match texture requirements");
    }

    return true;
}

void Material::loadImageFromFile(int index, const std::string& filename) {
    std::string fileExtension = filename.substr(filename.find_last_of(".") + 1);
    if (fileExtension == "png" || fileExtension == "PNG") {
        enabled[index] = loadImageFromPNG(texture[index], filename);
    } else if (fileExtension == "jpg" || fileExtension == "jpeg" || fileExtension == "JPG" || fileExtension == "JPEG") {
        enabled[index] = loadImageFromJPG(texture[index], filename);
    } else {
        //throw std::runtime_error("Unsupported file format");
    }
}