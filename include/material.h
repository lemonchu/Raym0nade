#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <string>
#include <assimp/material.h>
#include "geometry.h"

#define MAX_MIPMAP_LEVEL 8

class ImageData {
private:
    [[nodiscard]] int getPixelIndex(float u, float v, int level) const;
public:
    int width, height, channels;
    std::vector<uint8_t> data[MAX_MIPMAP_LEVEL];
    int map_depth;
    ImageData() : width(0), height(0), channels(0), map_depth(0) {}

    void generateMipmaps();
    [[nodiscard]] glm::vec3 get3(float u, float v, float depth) const;
    [[nodiscard]] glm::vec4 get4(float u, float v, float depth, bool gammaFlag) const;
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool hasTransparentPart() const;
};

std::string urlDecode(const std::string &src);

class Material {
public:

    ImageData texture[AI_TEXTURE_TYPE_MAX + 1];

    std::string name;
    bool hasTransparentPart;

    Material();

    [[nodiscard]] const ImageData& getImage(int index) const;

    static bool loadImageFromDDS(ImageData &imageData, const std::string& filename);
    static bool loadImageFromPNG(ImageData &imageData, const std::string& filename);
    static bool loadImageFromJPG(ImageData &imageData, const std::string& filename);
    void loadImageFromFile(int index, const std::string& filename);

    // New method to load material properties
    void loadMaterialProperties(const aiMaterial* aiMat);
    [[nodiscard]] vec4 getDiffuseColor(float u, float v) const;
    [[nodiscard]] vec3 getNormal(float u, float v) const;
    [[nodiscard]] vec3 getEmissiveColor(float u, float v) const;
    void getSurfaceData(float u, float v, float &rouguness, float &metallic) const;
};

#endif // MATERIAL_H