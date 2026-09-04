#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace raym0nade::detail {

struct ByteImage {
    int width{0};
    int height{0};
    int channels{0};
    std::vector<std::uint8_t> pixels;
};

struct FloatImage {
    int width{0};
    int height{0};
    int channels{0};
    std::vector<float> pixels;
};

[[nodiscard]] ByteImage loadByteImage(
    const std::filesystem::path& filename, const std::string& decoderModule);
[[nodiscard]] FloatImage loadFloatImage(
    const std::filesystem::path& filename, const std::string& decoderModule);

}  // namespace raym0nade::detail
