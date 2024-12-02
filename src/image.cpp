#include <vector>
#include <png.h>
#include <iostream>
#include "image.h"

RadianceData::RadianceData() {
    radiance = vec3(0.0f);
    Var = 0.0f;
    sampleCount = 0;
}

GbufferData::GbufferData() {
    diffuseColor = shapeNormal = surfaceNormal = position = emission = vec3(0.0f);
    face = nullptr;
    depth = 0.0f;
}

Image::Image(unsigned int width, unsigned int height) : width(width), height(height) {
    Gbuffer = new GbufferData[width * height];
    radiance_d = new RadianceData[width * height];
    radiance_i = new RadianceData[width * height];
    color = new vec3[width * height];
}

Image::~Image() {
    delete[] Gbuffer;
    delete[] radiance_d;
    delete[] radiance_i;
}

void filterVar(RadianceData *radiance, unsigned int width, unsigned int height) {
    static const int fliterRadius = 1;
    static const float h[3] = {0.25, 0.5, 0.25};
    std::vector<float> varBuffer(width * height, 0.0f);
    for (unsigned int y = 0; y < height; ++y)
        for (unsigned int x = 0; x < width; ++x)
            for (int dy = -fliterRadius; dy <= fliterRadius; ++dy)
                for (int dx = -fliterRadius; dx <= fliterRadius; ++dx) {
                    int nx = x + dx, ny = y + dy;
                    if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                        continue;
                    varBuffer[y * width + nx] += radiance[ny * width + nx].Var * h[dx + fliterRadius] * h[dy + fliterRadius];
                }
    for (unsigned int i = 0; i < width * height; ++i)
        radiance[i].Var = varBuffer[i];
}

float getWeight(const GbufferData &Gp, const GbufferData &Gq, const RadianceData &Lp, const RadianceData &Lq) {
    static const float
            sigma_z = 1.0f,
            sigma_n = 128.0f,
            sigma_l = 1.0f;

    if (Gq.face == nullptr)
        return 0.0f;

    float weight = 1.0f;

    vec3 dir = normalize(Gq.position - Gp.position);
    float sinTheta = abs(dot(Gp.shapeNormal, dir));
    float tanTheta = sinTheta / (sqrt(1.0f - sinTheta * sinTheta) + eps_zero);
    weight *= std::exp(- tanTheta / sigma_z);

    weight *= pow(std::max(0.0f, dot(Gp.surfaceNormal, Gq.surfaceNormal)), sigma_n);

    float radianceDiff = length(Lp.radiance - Lq.radiance) /
                           (sigma_l * sqrt(Lp.Var) + eps_zero);
    weight *= std::exp(-radianceDiff);

    return weight;
}

void filterRadiance(RadianceData *radiance, const GbufferData *Gbuffer, unsigned int width, unsigned int height, int step) {

    static const int filterRadius = 2;
    static const float h[5] = {0.0625, 0.25, 0.375, 0.25, 0.0625};
    std::vector<vec3> radianceBuffer(width * height, vec3(0.0f, 0.0f, 0.0f));
    std::vector<float> varBuffer(width * height, 0.0f);
    for (unsigned int y = 0; y < height; ++y)
        for (unsigned int x = 0; x < width; ++x) {
            float weightSum = 0.0f;
            RadianceData &Lp = radiance[y * width + x];
            const GbufferData &Gp = Gbuffer[y * width + x];
            if (Gp.face == nullptr)
                continue;
            for (int dy = -filterRadius; dy <= filterRadius; ++dy)
                for (int dx = -filterRadius; dx <= filterRadius; ++dx) {
                    int nx = x + dx * step, ny = y + dy * step;
                    if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                        continue;
                    float weight = h[dx + filterRadius] * h[dy + filterRadius];
                    RadianceData &Lq = radiance[ny * width + nx];
                    const GbufferData &Gq = Gbuffer[ny * width + nx];
                    if (dx != 0 || dy != 0)
                        weight *= getWeight(Gp, Gq, Lp, Lq);
                    weightSum += weight;
                    radianceBuffer[y * width + x] += Lq.radiance * weight;
                    varBuffer[y * width + x] += Lq.Var * weight * weight;
                }
            radianceBuffer[y * width + x] /= weightSum;
            varBuffer[y * width + x] /= weightSum * weightSum;
        }
    for (unsigned int i = 0; i < width * height; ++i) {
        radiance[i].radiance = radianceBuffer[i];
        radiance[i].Var = varBuffer[i];
    }
}

void filterRadiance(RadianceData *radiance, const GbufferData *Gbuffer, unsigned int width, unsigned int height) {

    filterVar(radiance, width, height);
    filterRadiance(radiance, Gbuffer, width, height, 1);
    filterRadiance(radiance, Gbuffer, width, height, 2);
    filterRadiance(radiance, Gbuffer, width, height, 4);
    filterRadiance(radiance, Gbuffer, width, height, 8);
    filterRadiance(radiance, Gbuffer, width, height, 16);
}

void normalize(RadianceData *radiance, unsigned int width, unsigned int height) {
    for (unsigned int i = 0; i < width * height; ++i) {
        radiance[i].radiance /= (float) radiance[i].sampleCount;
        radiance[i].Var /= (float) radiance[i].sampleCount;
        radiance[i].Var -= dot(radiance[i].radiance, radiance[i].radiance);
    }
}

void Image::normalizeRadiance() {
    normalize(radiance_d, width, height);
    normalize(radiance_i, width, height);
}

void Image::filter() {
    filterRadiance(radiance_d, Gbuffer, width, height);
    filterRadiance(radiance_i, Gbuffer, width, height);
}

void Image::shade() {
    for (unsigned int i = 0; i < width * height; ++i) {
        GbufferData &G = Gbuffer[i];
        RadianceData &Rd = radiance_d[i], &Ri = radiance_i[i];
        color[i] = G.diffuseColor * (Rd.radiance + Ri.radiance) + G.emission;
    }
}

void Image::bloom() {
    static const int filterRadius = 2;
    static const float h[5] = {0.0625, 0.25, 0.375, 0.25, 0.0625};
    std::vector<vec3> bloomColor(width * height);

    static const float omega_bloom = 0.2f, threshold = 1.5f;
    for (int id = 0; id < width * height; ++id) {
        bloomColor[id] = glm::max(color[id] - vec3(threshold), vec3(0.0f));
        bloomColor[id] = glm::pow(bloomColor[id], vec3(omega_bloom));
    }

    for (int step = 1; step <= 16; step *= 2) {
        std::vector<vec3> tempBloomColor = bloomColor;
        for (int y = 0; y < height; ++y)
            for (int x = 0; x < width; ++x) {
                vec3 &col = bloomColor[y * width + x];
                col = vec3(0.0f);
                for (int ky = -filterRadius; ky <= filterRadius; ++ky)
                    for (int kx = -filterRadius; kx <= filterRadius; ++kx) {
                        int nx = x + kx * step, ny = y + ky * step;
                        if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                            continue;
                        col += tempBloomColor[ny * width + nx] * h[ky + 2] * h[kx + 2];
                    }
            }
        for (unsigned int id = 0; id < width * height; ++id)
            color[id] += bloomColor[id] / 8.0f;
    }
}

const float gamma = 2.2f;

void Image::gammaCorrection() {
    for (unsigned int i = 0; i < width * height; ++i) {
        color[i] = glm::min(glm::max(color[i], vec3(0.0f) ), vec3(1.0f));
        color[i] = glm::pow(color[i], vec3(1.0f / gamma));
    }
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
            for (int k = 0; k < 3; ++k)
                row[k] = static_cast<png_byte>(color[id][k] * 255);
        }

    for (unsigned int y = 0; y < height; ++y)
        png_write_row(png_ptr, &image_data[y * width * 3]);

    png_write_end(png_ptr, nullptr);

    fclose(fp);
    png_destroy_write_struct(&png_ptr, &info_ptr);
}