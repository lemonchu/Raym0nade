#include <exception>
#include <filesystem>
#include <iostream>

#include "raym0nade/image.hpp"

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: raym0nade_fxaa <input.png> <output.png>\n";
        return 2;
    }

    try {
        raym0nade::Film image{std::filesystem::path{argv[1]}};
        image.reverseGammaCorrect();
        image.applyFxaa();
        image.gammaCorrect();
        image.save(std::filesystem::path{argv[2]});
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "FXAA failed: " << exception.what() << '\n';
        return 1;
    }
}
