#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <string>
#include <assimp/material.h>
#include <glm/glm.hpp>
#include <FreeImage.h>

struct ImageData{
    int width, height;
    std::vector<uint8_t> data;
    ImageData() : width(0), height(0) {}
};

std::string urlDecode(const std::string &src);

class Material {
public:

    ImageData texture[AI_TEXTURE_TYPE_MAX + 1];

    std::string name;
    float shininess;
    float opacity;
    glm::vec3 diffuseColor;
    glm::vec3 specularColor;
    glm::vec3 ambientColor;
    bool isNameEnabled;
    bool isShininessEnabled;
    bool isOpacityEnabled;
    bool isDiffuseColorEnabled;
    bool isSpecularColorEnabled;
    bool isAmbientColorEnabled;

    Material();

    [[nodiscard]] const ImageData& getImage(int index) const;

    int isEnabled(int index) const;

    bool loadImageFromDDS(ImageData &imageData, const std::string& filename);
    bool loadImageFromPNG(ImageData &imageData, const std::string& filename);
    bool loadImageFromJPG(ImageData &imageData, const std::string& filename);
    void loadImageFromFile(int index, const std::string& filename);

    // New method to load material properties
    void loadMaterialProperties(const aiMaterial* aiMat);
};

#endif // MATERIAL_H