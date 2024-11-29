#include <vector>
#include <png.h>
#include <iostream>
#include "image.h"

PixelData::PixelData() {
    face = nullptr;
    sampleCount = 0;
    Epower2 = 0.0f;
    inradiance = vec3(0.0f);
}

Image::Image(unsigned int width, unsigned int height) : buffer(nullptr), width(width), height(height) {
    buffer = new PixelData[width * height];
}

const int VarRadius = 2;
const float baseRatio = 0.25f;

/*void Image::allocateSamples(int totalSamples) {

    float totalNearVariance = 0.0f;

    for (unsigned int y = 0; y < height; ++y) {
        for (unsigned int x = 0; x < width; ++x) {
            float mean = 0.0f, mean2 = 0.0f;
            int xL = std::max(0, static_cast<int>(x) - VarRadius),
                xR = std::min(static_cast<int>(width), static_cast<int>(x) + VarRadius + 1),
                yL = std::max(0, static_cast<int>(y) - VarRadius),
                yR = std::min(static_cast<int>(height), static_cast<int>(y) + VarRadius + 1);
            for (int i = yL; i < yR; ++i)
                for (int j = xL; j < xR; ++j) {
                    vec4 color = buffer[i * width + j].color;
                    color[0] /= color[3];
                    color[1] /= color[3];
                    color[2] /= color[3];
                    float power = length(color);
                    mean += power;
                    mean2 += power * power;
                }
            mean /= (xR - xL) * (yR - yL);
            mean2 /= (xR - xL) * (yR - yL);
            buffer[y * width + x].VarNear = mean2 - mean * mean + eps_zero;
            totalNearVariance += buffer[y * width + x].VarNear;
        }
    }

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<float> U(0.0f, 1.0f);
    float baseProb = baseRatio * totalSamples / (width * height);
    for (unsigned int y = 0; y < height; ++y) {
        for (unsigned int x = 0; x < width; ++x) {
            float sampleCount = totalSamples * (1.0f - baseRatio) * (buffer[y * width + x].VarNear / totalNearVariance) + baseProb;
            int intSampleCount = static_cast<int>(sampleCount);
            float fractionalPart = sampleCount - intSampleCount;
            if (U(gen) < fractionalPart)
                intSampleCount++;
            buffer[y * width + x].samples = intSampleCount;
        }
    }

    std::cout << "Samples allocated." << std::endl;
}*/


/*void Image::filterInradiance(int step) {
    const float h[5] = {0.0625, 0.25, 0.375, 0.25, 0.0625};

}*/

void Image::filterVar() {
    const int fliterRadius = 1;
    const float h[3] = {0.25, 0.5, 0.25};
    std::vector<float> varBuffer(width * height, 0.0f);
    for (unsigned int y = 0; y < height; ++y)
        for (unsigned int x = 0; x < width; ++x)
            for (int dy = -fliterRadius; dy <= fliterRadius; ++dy)
                for (int dx = -fliterRadius; dx <= fliterRadius; ++dx) {
                    int nx = x + dx, ny = y + dy;
                    if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                        continue;
                    varBuffer[y * width + nx] += buffer[ny * width + nx].Var * h[dx + fliterRadius] * h[dy + fliterRadius];
                }
    for (unsigned int i = 0; i < width * height; ++i)
        buffer[i].Var = varBuffer[i];
}

float getWeight(const PixelData &current, const PixelData &neighbor) {
    const float
            sigma_z = 1.0f,
            sigma_n = 64.0f,
            sigma_l = 2.0f;

    if (neighbor.face == nullptr)
        return 0.0f;

    float weight = 1.0f;

    vec3 dir = normalize(neighbor.position - current.position);
    float sinTheta = abs(dot(current.shapeNormal, dir));
    float tanTheta = sinTheta / (sqrt(1.0f - sinTheta * sinTheta) + eps_zero);
    weight *= std::exp(- tanTheta / sigma_z);

    weight *= pow(std::max(0.0f, dot(current.surfaceNormal, neighbor.surfaceNormal)), sigma_n);

    float inradianceDiff = length(current.inradiance - neighbor.inradiance) /
                           (sigma_l * sqrt(neighbor.Var) + eps_zero);
    weight *= std::exp(-inradianceDiff);

    return weight;
}

void Image::filterInradiance(int step) {

    const int fliterRadius = 2;
    const float h[5] = {0.0625, 0.25, 0.375, 0.25, 0.0625};
    std::vector<vec3> inradianceBuffer(width * height, vec3(0.0f, 0.0f, 0.0f));
    std::vector<float> varBuffer(width * height, 0.0f);
    for (unsigned int y = 0; y < height; ++y)
        for (unsigned int x = 0; x < width; ++x) {
            float weightSum = 0.0f;
            PixelData &current = buffer[y * width + x];
            if (current.face == nullptr)
                continue;
            for (int dy = -fliterRadius; dy <= fliterRadius; ++dy)
                for (int dx = -fliterRadius; dx <= fliterRadius; ++dx) {
                    int nx = x + dx * step, ny = y + dy * step;
                    if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                        continue;
                    float weight = h[dx + fliterRadius] * h[dy + fliterRadius];
                    PixelData &neighbor = buffer[ny * width + nx];
                    if (dx != 0 || dy != 0)
                        weight *= getWeight(current, neighbor);
                    weightSum += weight;
                    inradianceBuffer[y * width + x] += neighbor.inradiance * weight;
                    varBuffer[y * width + x] += neighbor.Var * weight * weight;
                }
            inradianceBuffer[y * width + x] /= weightSum;
            varBuffer[y * width + x] /= weightSum * weightSum;
        }
    for (unsigned int i = 0; i < width * height; ++i) {
        buffer[i].inradiance = inradianceBuffer[i];
        buffer[i].Var = varBuffer[i];
    }
}

void Image::filterInradiance() {

    for (unsigned int i = 0; i < width * height; ++i) {
        buffer[i].inradiance /= (float) buffer[i].sampleCount;
        buffer[i].Epower2 /= (float) buffer[i].sampleCount;
        buffer[i].Var = buffer[i].Epower2 - dot(buffer[i].inradiance, buffer[i].inradiance);
    }

    filterVar();
    filterInradiance(1);
    filterInradiance(2);
    filterInradiance(4);
    filterInradiance(8);
    filterInradiance(16);

    for (unsigned int i = 0; i < width * height; ++i)
        buffer[i].inradiance += buffer[i].emission;
}


void Image::save(const char *file_name) {

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);

    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        std::cerr << "Error during png creation" << std::endl;
        return;
    }

    FILE *fp;
    if ((fp = fopen(file_name, "wb")) == nullptr) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        std::cerr << "Could not open file for writing" << std::endl;
        return;
    }

    png_init_io(png_ptr, fp);

    png_set_IHDR(
            png_ptr,
            info_ptr,
            width,
            height,
            8,
            PNG_COLOR_TYPE_RGB,
            PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT,
            PNG_FILTER_TYPE_DEFAULT
    );

    png_write_info(png_ptr, info_ptr);

    std::vector<png_byte> image_data;
    image_data.resize(width * height * 3);
    for (unsigned int y = 0; y < height; ++y)
        for (unsigned int x = 0; x < width; ++x) {
            int id = y * width + x;
            png_byte *row = &image_data[id * 3];
            PixelData &pixel = buffer[id];
            vec3 color = pixel.diffuseColor * pixel.inradiance;
            /*pixel.Epower /= pixel.sampleCount;
            pixel.Epower2 /= pixel.sampleCount;
            float Var = pixel.Epower2 - pixel.Epower * pixel.Epower;
            row[0] = row[1] = row[2] = static_cast<png_byte>(Var * 128);*/
            for (int k = 0; k < 3; ++k) {
                float buf = std::min(std::max(color[k], 0.0f), 1.0f);
                buf = std::pow(buf,1.0f / gamma);
                row[k] = static_cast<png_byte>(buf * 255);
            }

        }

    for (unsigned int y = 0; y < height; ++y)
        png_write_row(png_ptr, &image_data[y * width * 3]);

    png_write_end(png_ptr, nullptr);

    fclose(fp);
    png_destroy_write_struct(&png_ptr, &info_ptr);
}

Image::~Image() {
    delete[] buffer;
}