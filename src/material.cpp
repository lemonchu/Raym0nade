#include <assimp/material.h>
#include <cstdlib>
#include <png.h>
#include <turbojpeg.h>
#include <stdexcept>
#include <vector>
#include <iostream>
#include "material.h"

#include <sstream>

std::string urlDecode(const std::string &src) {
    std::string decoded;
    char ch;
    int i, ii;
    for (i = 0; i < src.length(); i++) {
        if (int(src[i]) == 37) {
            sscanf(src.substr(i + 1, 2).c_str(), "%x", &ii);
            ch = static_cast<char>(ii);
            decoded += ch;
            i = i + 2;
        } else {
            decoded += src[i];
        }
    }
    return decoded;
}

glm::vec4 ImageData::get(float u, float v) const {
    int texX = lround(u * width), texY = lround(v * height);
    texX %= width;
    if (texX < 0) texX += width;
    texY %= height;
    if (texY < 0) texY += height;
    int pixelIndex = (texY * width + texX) * 4; // 4 channels for RGBA
    return glm::vec4(
            gammaMap[data[pixelIndex]],
            gammaMap[data[pixelIndex + 1]],
            gammaMap[data[pixelIndex + 2]],
            data[pixelIndex + 3] / 255.0f
    );
}

bool ImageData::empty() const {
    return data.empty();
}

bool ImageData::hasTransparentPart() const {
    for (int i = 3; i < data.size(); i += 4)
        if (data[i] < 255)
            return true;
    return false;
}

Material::Material() : shininess(0.0f), opacity(1.0f),
                       isNameEnabled(false), isShininessEnabled(false), isOpacityEnabled(false),
                       isDiffuseColorEnabled(false), isSpecularColorEnabled(false), isAmbientColorEnabled(false),
                       isEmissionEnabled(false) {
}

[[nodiscard]] const ImageData& Material::getImage(int index) const {
    if (index < 0 || index >= AI_TEXTURE_TYPE_MAX + 1) {
        throw std::out_of_range("Index out of range");
    }
    return texture[index];
}

bool Material::loadImageFromDDS(ImageData &imageData, const std::string& filename) {
    FIBITMAP* bitmap = FreeImage_Load(FIF_DDS, filename.c_str(), DDS_DEFAULT);
    if (!bitmap) {
        std::cerr << "Error loading DDS file: " << filename << std::endl;
        return false;
    }

    unsigned int imgWidth = FreeImage_GetWidth(bitmap);
    unsigned int imgHeight = FreeImage_GetHeight(bitmap);
    std::cout << "Loading DDS image: " << imgWidth << "x" << imgHeight << std::endl;
    unsigned int pitch = imgWidth * 4; // Assuming 4 bytes per pixel (RGBA)
    imageData.data.resize(imgHeight * pitch);

    for (int y = 0; y < imgHeight; y++) {
        BYTE* bits = FreeImage_GetScanLine(bitmap, y);
        for (int x = 0; x < imgWidth; x++) {
            imageData.data[y * pitch + 4 * x + 0] = bits[4 * x + 2]; // Red
            imageData.data[y * pitch + 4 * x + 1] = bits[4 * x + 1]; // Green
            imageData.data[y * pitch + 4 * x + 2] = bits[4 * x + 0]; // Blue
            imageData.data[y * pitch + 4 * x + 3] = bits[4 * x + 3]; // Alpha
        }
    }
    imageData.width = imgWidth;
    imageData.height = imgHeight;

    FreeImage_Unload(bitmap);
    return true;
}

bool Material::loadImageFromPNG(ImageData &imageData, const std::string& filename) {

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
    imageData.data.resize(imgWidth * imgHeight * 3);

    for (int y = 0; y < imgHeight; y++) {
        png_read_row(png, row.data(), nullptr);
        for (int x = 0; x < imgWidth * 3; x++) {
            imageData.data[y * imgWidth * 3 + x] = row[x];
        }
    }
    imageData.width = imgWidth;
    imageData.height = imgHeight;

    fclose(fp);
    png_destroy_read_struct(&png, &info, nullptr);

    return true;
}

bool Material::loadImageFromJPG(ImageData &imageData, const std::string& filename) {
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
    imageData.data.resize(imgHeight * pitch);

    // Decompress JPEG image to RGB format
    if (tjDecompress2(tjInstance, jpegBuf, jpegSize, &imageData.data[0], imgWidth, pitch, imgHeight, TJPF_RGB, TJFLAG_FASTDCT) < 0) {
        std::cerr << "Error decompressing JPEG image" << std::endl;
        free(jpegBuf);
        tjDestroy(tjInstance);
        return false;
    }
    imageData.width = imgWidth;
    imageData.height = imgHeight;

    // Free resources
    free(jpegBuf);
    tjDestroy(tjInstance);

    return true;
}

void Material::loadImageFromFile(int index, const std::string& filename) {
    std::string fileExtension = filename.substr(filename.find_last_of(".") + 1);
    if (fileExtension == "png" || fileExtension == "PNG") {
        loadImageFromPNG(texture[index], filename);
    } else if (fileExtension == "jpg" || fileExtension == "jpeg" || fileExtension == "JPG" || fileExtension == "JPEG") {
        loadImageFromJPG(texture[index], filename);
    } else if (fileExtension == "dds" || fileExtension == "DDS") {
        loadImageFromDDS(texture[index], filename);
    } else {
        std::cerr << "Unsupported file format: " << filename << std::endl;
    }
}

void Material::loadMaterialProperties(const aiMaterial* aiMat) {
    aiString matName;
    if (aiMat->Get(AI_MATKEY_NAME, matName) == AI_SUCCESS) {
        name = matName.C_Str();
        isNameEnabled = true;
    }

    if (aiMat->Get(AI_MATKEY_SHININESS, shininess) == AI_SUCCESS) {
        isShininessEnabled = true;
    }

    if (aiMat->Get(AI_MATKEY_OPACITY, opacity) == AI_SUCCESS) {
        isOpacityEnabled = true;
    }

    aiColor3D color;
    if (aiMat->Get(AI_MATKEY_COLOR_DIFFUSE, color) == AI_SUCCESS) {
        diffuseColor = glm::vec3(color.r, color.g, color.b);
        isDiffuseColorEnabled = true;
    }

    if (aiMat->Get(AI_MATKEY_COLOR_SPECULAR, color) == AI_SUCCESS) {
        specularColor = glm::vec3(color.r, color.g, color.b);
        isSpecularColorEnabled = true;
    }

    if (aiMat->Get(AI_MATKEY_COLOR_AMBIENT, color) == AI_SUCCESS) {
        ambientColor = glm::vec3(color.r, color.g, color.b);
        isAmbientColorEnabled = true;
    }

    if (aiMat->Get(AI_MATKEY_COLOR_EMISSIVE, color) == AI_SUCCESS) {
        emission = glm::vec3(color.r, color.g, color.b);
        if (color.r > 0.0f || color.g > 0.0f || color.b > 0.0f)
            isEmissionEnabled = true;
    }

    if (texture[TextureIdForDiffuseColor].hasTransparentPart()) {
        hasTransparentPart = true;
        std::cout << "Material has transparent part" << std::endl;
    }
}

glm::vec4 Material::getDiffuseColor(float u, float v) const {
    if (texture[TextureIdForDiffuseColor].empty()) {
        return glm::vec4(diffuseColor, 1.0f);
    }
    return texture[TextureIdForDiffuseColor].get(u, v);
}