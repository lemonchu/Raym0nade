#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <string>
#include <assimp/material.h>
#include "geometry.h"

#define MAX_MIPMAP_LEVEL 8

class ImageData {
public:
    int width, height, channels;
    std::vector<uint8_t> data[MAX_MIPMAP_LEVEL];
    int map_depth;
    ImageData();

    void generateMipmaps();
    template<typename Vec>
    [[nodiscard]] Vec get(float u, float v, float depth) const;
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool hasTransparentPart() const;
};

class ImageDataHDR {
private:
    void InitSunDir();
public:
    int width, height;
    std::vector<float> data;
    vec3 sunDir;
    ImageDataHDR();
    void load(const std::string &filename);
    [[nodiscard]] bool empty() const;
    [[nodiscard]] vec3 get(vec3 dir) const;
};

std::string urlDecode(const std::string &src);

class Material {
public:

    ImageData texture[AI_TEXTURE_TYPE_MAX];

    std::string name;
    int id;
    bool hasFullyTransparentPart;
    float opacity, ior, roughness;
    vec3 transmittingColor;

    Material();

    static bool loadImageFromDDS(ImageData &imageData, const std::string& filename);
    static bool loadImageFromPNG(ImageData &imageData, const std::string& filename);
    void loadImageFromFile(int index, const std::string& filename);

    // New method to load material properties
    void loadMaterialProperties(const aiMaterial* aiMat);
    [[nodiscard]] vec4 getDiffuseColor(float u, float v, float duv) const;
    [[nodiscard]] vec3 getNormal(float u, float v, float duv) const;
    [[nodiscard]] vec3 getEmissiveColor(float u, float v, float duv) const;
    void getSurfaceData(float u, float v, float &roughness, float &metallic) const;
};

#endif // MATERIAL_H