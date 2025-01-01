#include <cstdlib>
#include <png.h>
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

template<typename Vec>
Vec get_basic(const std::vector<uint8_t> &data, int index);

template<>
vec3 get_basic<vec3>(const std::vector<uint8_t> &data, int index) {
    return vec3(
            static_cast<float>(data[index]),
            static_cast<float>(data[index + 1]),
            static_cast<float>(data[index + 2])
    ) / 255.0f;
}

template<>
vec4 get_basic<vec4>(const std::vector<uint8_t> &data, int index) {
    return vec4(
            static_cast<float>(data[index]),
            static_cast<float>(data[index + 1]),
            static_cast<float>(data[index + 2]),
            static_cast<float>(data[index + 3])
    ) / 255.0f;
}

template<typename Vec>
Vec get_bilinear(const std::vector<uint8_t>& data, int width, int height, float u, float v) {

    static constexpr auto _mod = [](int &x, int m) -> void {
        if (x < 0 || x >= m) {
            x %= m;
            if (x < 0) x += m;
        }
    };

    static const int BitWidth = std::is_same<Vec, vec3>::value ? 3 : 4;
    float x = u * float(width), y = v * float(height);
    int x0 = int(floorf(x)), y0 = int(floorf(y));
    float dx = x - float(x0), dy = y - float(y0);
    _mod(x0, width);
    _mod(y0, height);
    int x1 = (x0 + 1) % width, y1 = (y0 + 1) % height;

    Vec c00 = get_basic<Vec>(data, (y0 * width + x0) * BitWidth);
    Vec c01 = get_basic<Vec>(data, (y0 * width + x1) * BitWidth);
    Vec c10 = get_basic<Vec>(data, (y1 * width + x0) * BitWidth);
    Vec c11 = get_basic<Vec>(data, (y1 * width + x1) * BitWidth);
    Vec c0 = c00 * (1.0f - dx) + c01 * dx;
    Vec c1 = c10 * (1.0f - dx) + c11 * dx;

    Vec ret = c0 * (1.0f - dy) + c1 * dy;

    if (dx < 0.0f || dy < 0.0f)
        std::cout << "dx: " << dx << " dy: " << dy << std::endl;

    return c0 * (1.0f - dy) + c1 * dy;
}

template<typename Vec>
Vec ImageData::get(float u, float v, float depth) const {
    v = 1.0f - v;
    depth = std::max(std::min(depth, static_cast<float>(map_depth - 1)), 0.0f);

    int level = static_cast<int>(depth);
    int nextLevel = std::min(level + 1, map_depth - 1);
    float levelBlend = depth - float(level);

    Vec ret1 = get_bilinear<Vec>(data[level], width >> level, height >> level, u, v);
    Vec ret2 = get_bilinear<Vec>(data[nextLevel], width >> nextLevel, height >> nextLevel, u, v);

    return ret1 * (1.0f - levelBlend) + ret2 * levelBlend;
}

ImageData::ImageData() : width(0), height(0), channels(0), map_depth(0) {}

bool ImageData::empty() const {
    return data[0].empty();
}

bool ImageData::hasTransparentPart() const {
    for (int i = 3; i < data[0].size(); i += 4)
        if (data[0][i] < 255)
            return true;
    return false;
}

Material::Material() :
    hasFullyTransparentPart(false), opacity(1.0f), ior(1.0f),
    transmittingColor(0.0f), roughness(0.8f), id(0) {}

void ImageData::generateMipmaps() {

    map_depth = MAX_MIPMAP_LEVEL; // NEVER REMOVE THIS LINE

    std::cout << "--- Generating mipmaps...";
    for (int level = 1; level < MAX_MIPMAP_LEVEL; ++level) {
        int prevWidth = width >> (level - 1);
        int prevHeight = height >> (level - 1);
        int currWidth = width >> level;
        int currHeight = height >> level;

        if (currWidth == 0 || currHeight == 0){
            map_depth = level;
            break;
        }

        data[level].resize(currWidth * currHeight * channels);

        for (int y = 0; y < currHeight; ++y) {
            for (int x = 0; x < currWidth; ++x) {
                for (int c = 0; c < channels; ++c) {
                    int sum = 0;
                    for (int dy = 0; dy < 2; ++dy) {
                        for (int dx = 0; dx < 2; ++dx) {
                            int px = (x * 2 + dx) % prevWidth;
                            int py = (y * 2 + dy) % prevHeight;
                            sum += data[level - 1][(py * prevWidth + px) * channels + c];
                        }
                    }
                    data[level][(y * currWidth + x) * channels + c] = sum / 4;
                }
            }
        }
    }
    std::cout << map_depth << " levels generated." << std::endl;
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
                uint8_t *buffer = static_cast<uint8_t*>(view.buf);
                imageData.data[0].assign(buffer, buffer + view.len);
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
    imageData.data[0].resize(imgWidth * imgHeight * 3);

    for (int y = 0; y < imgHeight; y++) {
        png_read_row(png, row.data(), nullptr);
        for (int x = 0; x < imgWidth * 3; x++) {
            imageData.data[0][y * imgWidth * 3 + x] = row[x];
        }
    }
    imageData.width = imgWidth;
    imageData.height = imgHeight;

    fclose(fp);
    png_destroy_read_struct(&png, &info, nullptr);

    return true;
}

void Material::loadImageFromFile(int index, const std::string &filename) {
    std::string fileExtension = filename.substr(filename.find_last_of(".") + 1);
    if (fileExtension == "png" || fileExtension == "PNG") {
        loadImageFromPNG(texture[index], filename);
    } else if (fileExtension == "dds" || fileExtension == "DDS") {
        loadImageFromDDS(texture[index], filename);
    } else {
        std::cerr << "Unsupported file format: " << filename << std::endl;
    }
    texture[index].generateMipmaps();
}

void Material::loadMaterialProperties(const aiMaterial *aiMat) {
    aiString matName;
    if (aiMat->Get(AI_MATKEY_NAME, matName) == AI_SUCCESS) {
        name = matName.C_Str();
    }

    ai_real aiOpacity;
    if (aiMat->Get(AI_MATKEY_OPACITY, aiOpacity) == AI_SUCCESS && aiOpacity < 0.99f) {
        opacity = 0.0f;
        ior = 1.25f;
        std::cout << "Material has Opacity" << std::endl;
        std::cout << "Opacity " << opacity << " ior: " << ior << std::endl;

        roughness = 5e-3;
        if (name == "Ice")
            roughness = 0.5f;

        if (name == "TransparentGlassWine")
            transmittingColor = vec3(0.2f, 0.08f, 0.07f);
        if (name == "TransparentGlass" || name == "Water" || name == "Ice")
            transmittingColor = vec3(1.0f);

        if (name == "Beer")
            transmittingColor = vec3(0.8f, 0.7f, 0.55f);
        if (name == "Red_Wine")
            transmittingColor = vec3(0.24f, 0.09f, 0.07f);
        if (name == "White_Wine")
            transmittingColor = vec3(0.85f, 0.78f, 0.6f);
    }

    if (texture[aiTextureType_DIFFUSE].hasTransparentPart()) {
        hasFullyTransparentPart = true;
        std::cout << "Material has transparent part" << std::endl;
    }

}

void gammaPow(vec4 &v) {
    static const float gamma = 2.2f;
    if (v[0] < 0.0f || v[1] < 0.0f || v[2] < 0.0f) {
        std::cerr << "color: " << v[0] << ", " << v[1] << ", " << v[2] << std::endl;
        throw std::runtime_error("Invalid color value");
    }
    v[0] = pow(v[0], gamma);
    v[1] = pow(v[1], gamma);
    v[2] = pow(v[2], gamma);
}

// 传入 duv = nan 表示不使用mipmap
vec4 Material::getDiffuseColor(float u, float v, float duv) const {
    if (texture[aiTextureType_DIFFUSE].empty())
        return vec4(1.0f);
    float map_depth = isnan(duv) ? 0.0f : log2(duv * float(texture[aiTextureType_DIFFUSE].width));
    vec4 color = texture[aiTextureType_DIFFUSE].get<vec4>(u, v, map_depth);
    gammaPow(color);
    return color;
}

vec3 Material::getNormal(float u, float v, float duv) const {
    if (texture[aiTextureType_NORMALS].empty())
        return vec3(0.0f);
    float map_depth = isnan(duv) ? 0.0f : log2(duv * float(texture[aiTextureType_NORMALS].width));
    return texture[aiTextureType_NORMALS].get<vec3>(u, v, map_depth) * 2.0f - 1.0f;
}

vec3 Material::getEmissiveColor(float u, float v, float duv) const {
    if (texture[aiTextureType_EMISSIVE].empty())
        return vec3(0.0f);
    float map_depth = isnan(duv) ? 0.0f : log2(duv * float(texture[aiTextureType_EMISSIVE].width));
    vec4 color = texture[aiTextureType_EMISSIVE].get<vec4>(u, v, map_depth);
    gammaPow(color);
    return color;
}

void Material::getSurfaceData(float u, float v, float &roughness, float &metallic) const {
    if (texture[aiTextureType_SPECULAR].empty()) {
        metallic = 0.0f;
        roughness = this->roughness;
        return;
    }
    vec4 surfaceData = texture[aiTextureType_SPECULAR].get<vec4>(u, v, 0);
    roughness = std::max(surfaceData[1], 1e-3f);
    metallic = std::min(surfaceData[2], 0.99f);
    roughness *= 0.5f + 0.5f * (1.0f - metallic);
}