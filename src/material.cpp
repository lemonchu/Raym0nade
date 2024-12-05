#include <cstdlib>
#include <png.h>
#include <turbojpeg.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include "material.h"
#include <Python.h>

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

int ImageData::getPixelIndex(float u, float v, int bitDepth) const {
    v = 1.0f - v;
    int texX = lround(u * width), texY = lround(v * height);
    if (texX < 0 || texX >= width) {
        texX %= width;
        if (texX < 0) texX += width;
    }
    if (texY < 0 || texY >= height) {
        texY %= height;
        if (texY < 0) texY += height;
    }
    return (texY * width + texX) * bitDepth;
}

glm::vec3 ImageData::get3(float u, float v) const {
    int pixelIndex = getPixelIndex(u, v, 3);
    return {
            static_cast<float>(data[pixelIndex]) / 255.0f,
            static_cast<float>(data[pixelIndex + 1]) / 255.0f,
            static_cast<float>(data[pixelIndex + 2]) / 255.0f
    };
}

glm::vec4 ImageData::get4(float u, float v) const {
    int pixelIndex = getPixelIndex(u, v, 4);
    return {
            static_cast<float>(data[pixelIndex]) / 255.0f,
            static_cast<float>(data[pixelIndex + 1]) / 255.0f,
            static_cast<float>(data[pixelIndex + 2]) / 255.0f,
            static_cast<float>(data[pixelIndex + 3]) / 255.0f
    };
}

glm::vec4 ImageData::get4_with_gamma(float u, float v) const {
    int pixelIndex = getPixelIndex(u, v, 4); // 4 channels for RGBA
    return {
            gammaMap[data[pixelIndex]],
            gammaMap[data[pixelIndex + 1]],
            gammaMap[data[pixelIndex + 2]],
            static_cast<float>(data[pixelIndex + 3]) / 255.0f
    };
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
                       isEmissionEnabled(false) {}

[[nodiscard]] const ImageData &Material::getImage(int index) const {
    if (index < 0 || index >= AI_TEXTURE_TYPE_MAX + 1) {
        throw std::out_of_range("Index out of range");
    }
    return texture[index];
}

PyObject *call_dds_to_array(const std::string &filename) {
    static bool isInitialized = false;

    // Initialize the Python interpreter only once
    if (!isInitialized) {
        Py_Initialize();
        isInitialized = true;

        // Add the scripts directory to the Python path
        PyObject *sysPath = PySys_GetObject("path");
        PyObject *path = PyUnicode_DecodeFSDefault("scripts");
        PyList_Append(sysPath, path);
        Py_DECREF(path);
    }

    PyObject *pName = PyUnicode_DecodeFSDefault("dds_to_array");
    PyObject *pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != nullptr) {
        PyObject *pFunc = PyObject_GetAttrString(pModule, "dds_to_array");
        if (PyCallable_Check(pFunc)) {
            PyObject *pArgs = PyTuple_Pack(1, PyUnicode_FromString(filename.c_str()));
            PyObject *pValue = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            Py_XDECREF(pFunc);
            Py_DECREF(pModule);
            return pValue;
        } else {
            PyErr_Print();
            Py_XDECREF(pFunc);
            Py_DECREF(pModule);
        }
    } else {
        PyErr_Print();
    }

    return nullptr;
}

bool Material::loadImageFromDDS(ImageData &imageData, const std::string &filename) {
    PyObject *pValue = call_dds_to_array(filename);

    if (pValue != nullptr) {
        PyObject *pWidth = PyTuple_GetItem(pValue, 0);
        PyObject *pHeight = PyTuple_GetItem(pValue, 1);
        PyObject *pChannels = PyTuple_GetItem(pValue, 2);
        PyObject *pData = PyTuple_GetItem(pValue, 3);

        if (pWidth != Py_None && pHeight != Py_None && pChannels != Py_None && pData != Py_None) {
            imageData.width = PyLong_AsLong(pWidth);
            imageData.height = PyLong_AsLong(pHeight);
            imageData.channels = PyLong_AsLong(pChannels);

            Py_buffer view;
            if (PyObject_GetBuffer(pData, &view, PyBUF_SIMPLE) == 0) {
                uint8_t *buffer = static_cast<uint8_t *>(view.buf);
                imageData.data.assign(buffer, buffer + view.len);
                PyBuffer_Release(&view);
            }
        }
        Py_DECREF(pValue);
    } else {
        PyErr_Print();
        return false;
    }

    return true;
}

bool Material::loadImageFromPNG(ImageData &imageData, const std::string &filename) {

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

bool Material::loadImageFromJPG(ImageData &imageData, const std::string &filename) {
    // Open JPEG file
    FILE *jpegFile = fopen(filename.c_str(), "rb");
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
    unsigned char *jpegBuf = (unsigned char *) malloc(jpegSize);
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
    if (tjDecompress2(tjInstance, jpegBuf, jpegSize, &imageData.data[0], imgWidth, pitch, imgHeight, TJPF_RGB,
                      TJFLAG_FASTDCT) < 0) {
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

void Material::loadImageFromFile(int index, const std::string &filename) {
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

void Material::loadMaterialProperties(const aiMaterial *aiMat) {
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

    if (texture[aiTextureType_DIFFUSE].hasTransparentPart()) {
        hasTransparentPart = true;
        std::cout << "Material has transparent part" << std::endl;
    }
}

glm::vec4 Material::getDiffuseColor(float u, float v) const {
    if (texture[aiTextureType_DIFFUSE].empty())
        return vec4(diffuseColor, 1.0f);
    return texture[aiTextureType_DIFFUSE].get4_with_gamma(u, v);
}

glm::vec3 Material::getNormal(float u, float v) const {
    if (texture[aiTextureType_NORMALS].empty())
        return vec3(0.0f, 0.0f, 1.0f);
    return 1.0f - texture[aiTextureType_NORMALS].get3(u, v) * 2.0f;
}