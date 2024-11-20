#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <string>
#include <assimp/material.h>
#include <glm/glm.hpp>

class Material {
public:
    int width, height;
    std::vector<uint8_t> texture[AI_TEXTURE_TYPE_MAX];
    bool enabled[AI_TEXTURE_TYPE_MAX]; // Array to indicate if each texture is enabled

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

    [[nodiscard]] const std::vector<uint8_t>& getImage(int index) const;

    bool loadImageFromPNG(std::vector<uint8_t> &imageData, const std::string& filename);
    bool loadImageFromJPG(std::vector<uint8_t> &imageData, const std::string& filename);
    void loadImageFromFile(int index, const std::string& filename);

    // New method to load material properties
    void loadMaterialProperties(const aiMaterial* aiMat);
};

#endif // MATERIAL_H