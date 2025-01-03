#include <iostream>
#include "image.h"
#include "glm/glm.hpp"

int main() {
    const char* file_name = "./output/raw_2.png";
    Image png_to_be_fxaa(file_name);

    png_to_be_fxaa.reverseGammaCorrection();
    png_to_be_fxaa.FXAA();
    png_to_be_fxaa.gammaCorrection();

    png_to_be_fxaa.save("./output/raw_FXAA.png");

    return 0;
}