#ifndef TEXTURE_H
#define TEXTURE_H

#include <vector>
#include <assimp/material.h>

class Texture {
public:
    int width, height;
    std::vector<uint8_t> images[AI_TEXTURE_TYPE_MAX];
    bool enabled[AI_TEXTURE_TYPE_MAX]; // Array to indicate if each image is enabled

    Texture();

    [[nodiscard]] const std::vector<uint8_t>& getImage(int index) const;

    bool loadImageFromPNG(std::vector<uint8_t> &imageData, const std::string& filename);
    bool loadImageFromJPG(std::vector<uint8_t> &imageData, const std::string& filename);
    void loadImageFromFile(int index, const std::string& filename);
};

#endif // TEXTURE_H