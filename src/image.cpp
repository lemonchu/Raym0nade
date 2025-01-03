#include <vector>
#include <png.h>
#include <iostream>
#include <thread>
#include "image.h"

RadianceData::RadianceData() : radiance(vec3(0.0f)), Var(0.0f) {}

Image::Image(int width, int height) : width(width), height(height) {
    Gbuffer = new HitInfo[width * height];
    radiance_Dd = new RadianceData[width * height];
    radiance_Ds = new RadianceData[width * height];
    radiance_Id = new RadianceData[width * height];
    radiance_Is = new RadianceData[width * height];
    pixelarray = new vec3[width * height];
}

Image::Image(const char *file_name): radiance_Dd(nullptr), radiance_Ds(nullptr), radiance_Id(nullptr), radiance_Is(nullptr), pixelarray(nullptr), Gbuffer(nullptr) {
    load(file_name);
}

Image::~Image() {
    delete[] Gbuffer;
    delete[] radiance_Dd;
    delete[] radiance_Ds;
    delete[] radiance_Id;
    delete[] radiance_Is;
}

void filterVar(RadianceData *radiance, int width, int height) {
    static const int fliterRadius = 1;
    static const float h[3] = {0.25, 0.5, 0.25};
    std::vector<float> varBuffer(width * height, 0.0f);
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            vec3 E = vec3(0.0f);
            float E2 = 0.0f, Var = 0.0f;
            for (int dy = -fliterRadius; dy <= fliterRadius; ++dy)
                for (int dx = -fliterRadius; dx <= fliterRadius; ++dx) {
                    int nx = x + dx, ny = y + dy;
                    if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                        continue;
                    RadianceData &Lq = radiance[ny * width + nx];
                    float weight = h[dx + fliterRadius] * h[dy + fliterRadius];
                    E += Lq.radiance * weight;
                    E2 += dot(Lq.radiance, Lq.radiance) * weight;
                    Var += Lq.Var * weight;
                }
            varBuffer[y * width + x] = Var + E2 - dot(E, E);
        }
    for (int i = 0; i < width * height; ++i)
        radiance[i].Var = varBuffer[i];
}

float getWeight(const HitInfo &Gp, const HitInfo &Gq, const RadianceData &Lp, const RadianceData &Lq) {
    static constexpr float
            sigma_z = 1.0f,
            sigma_n = 1024.0f,
            sigma_l = 1.0f,
            sigma_m = 1.0f,
            eps = 1e-2f;

    if (!isfinite(Gq.position))
        return 0.0f;

    float weight = pow(std::max(0.0f, dot(Gp.surfaceNormal, Gq.surfaceNormal)), sigma_n);

    if (weight < 1e-6f)
        return 0.0f;

    float k = 0.0f;

    vec3 dir = normalize(Gq.position - Gp.position);
    float sinTheta = abs(dot(Gp.shapeNormal, dir));
    float tanTheta = sinTheta / (sqrt(1.0f - sinTheta * sinTheta) + eps_zero);
    k += - tanTheta / sigma_z;

    float radianceDiff = length(Lp.radiance - Lq.radiance) /
                           (sigma_l * sqrt(Lp.Var) + eps);
    k += - radianceDiff;

    float materialDiff = length(vec3(
        Gp.metallic - Gq.metallic,
        Gp.specular - Gq.specular,
        Gp.opacity - Gq.opacity
    ));
    k += - materialDiff / sigma_m;

    if (k < -7.5f)
        return 0.0f;

    weight *= exp(k);

    weight *= std::max(Gp.roughness, 0.04f);

    if (!isfinite(weight))
        return 0.0f;

    return weight;
}

void filterRadiance(RadianceData *radiance, const HitInfo *Gbuffer, int width, int height, int step) {

    static const int filterRadius = 2;
    static const float h[5] = {0.0625, 0.25, 0.375, 0.25, 0.0625};
    std::vector<vec3> radianceBuffer(width * height, vec3(0.0f, 0.0f, 0.0f));
    std::vector<float> varBuffer(width * height, 0.0f);
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            float weightSum = 0.0f;
            RadianceData &Lp = radiance[y * width + x];
            const HitInfo &Gp = Gbuffer[y * width + x];
            if (!isfinite(Gp.position))
                continue;
            for (int dy = -filterRadius; dy <= filterRadius; ++dy)
                for (int dx = -filterRadius; dx <= filterRadius; ++dx) {
                    int nx = x + dx * step, ny = y + dy * step;
                    if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                        continue;
                    float weight = h[dx + filterRadius] * h[dy + filterRadius];
                    RadianceData &Lq = radiance[ny * width + nx];
                    const HitInfo &Gq = Gbuffer[ny * width + nx];
                    if (dx != 0 || dy != 0)
                        weight *= getWeight(Gp, Gq, Lp, Lq);
                    weightSum += weight;
                    radianceBuffer[y * width + x] += Lq.radiance * weight;
                    varBuffer[y * width + x] += Lq.Var * weight * weight;
                }
            radianceBuffer[y * width + x] /= weightSum;
            varBuffer[y * width + x] /= weightSum * weightSum;
        }
    for (int i = 0; i < width * height; ++i) {
        radiance[i].radiance = radianceBuffer[i];
        radiance[i].Var = varBuffer[i];
    }
}

void filterRadiance(RadianceData *radiance, const HitInfo *Gbuffer, int width, int height) {
    filterVar(radiance, width, height);
    filterRadiance(radiance, Gbuffer, width, height, 1);
    filterRadiance(radiance, Gbuffer, width, height, 2);
    filterRadiance(radiance, Gbuffer, width, height, 4);
    filterRadiance(radiance, Gbuffer, width, height, 8);
    filterRadiance(radiance, Gbuffer, width, height, 16);
}

void Image::filter() {
    std::thread t1([this]() { filterRadiance(radiance_Dd, Gbuffer, width, height); });
    std::thread t2([this]() { filterRadiance(radiance_Ds, Gbuffer, width, height); });
    std::thread t3([this]() { filterRadiance(radiance_Id, Gbuffer, width, height); });
    std::thread t4([this]() { filterRadiance(radiance_Is, Gbuffer, width, height); });

    t1.join();
    t2.join();
    t3.join();
    t4.join();
}

void Image::shade(float exposure, int options) {
    for (int i = 0; i < width * height; ++i) {
        const HitInfo &G = Gbuffer[i];
        if (options & shapeNormal) {
            pixelarray[i] = (G.shapeNormal + vec3(1.0f)) / 2.0f;
            continue;
        }
        if (options & surfaceNormal) {
            pixelarray[i] = (G.surfaceNormal + vec3(1.0f)) / 2.0f;
            continue;
        }
        vec3 radiance_d = vec3(0.0f), radiance_s = vec3(0.0f);
        if (options & DirectLight) {
            if (options & Diffuse)
                radiance_d += radiance_Dd[i].radiance;
            if (options & Specular)
                radiance_s += radiance_Ds[i].radiance;
        }
        if (options & IndirectLight) {
            if (options & Diffuse)
                radiance_d += radiance_Id[i].radiance;
            if (options & Specular)
                radiance_s += radiance_Is[i].radiance;
        }
        if (!(options & (DirectLight | IndirectLight)))
            radiance_d = vec3(1.0f);
        vec3 diffuseColor = (options & BaseColor) ? G.baseColor : vec3(1.0f);
        pixelarray[i] = diffuseColor * radiance_d + radiance_s;
        if (options & Emission)
            pixelarray[i] += G.emission * exposure;
    }
}

void Image::bloom() {
    static const int filterRadius = 2;
    static const float h[5] = {0.0625, 0.25, 0.375, 0.25, 0.0625};
    std::vector<vec3> bloomColor(width * height);

    static const float omega_bloom = 0.85f;
    for (int id = 0; id < width * height; ++id) {
        float Clum = dot(pixelarray[id], RGB_Weight);
        if (Clum < 1.0f)
            continue;
        bloomColor[id] = pixelarray[id] / pow(Clum, omega_bloom);
        for (int k = 0; k < 3; ++k)
            bloomColor[id][k] = std::max(bloomColor[id][k]-1.0f, 0.0f);
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
                        col += tempBloomColor[ny * width + nx] * h[ky + filterRadius] * h[kx + filterRadius];
                    }
            }
        for (int id = 0; id < width * height; ++id)
            pixelarray[id] += bloomColor[id] / 8.0f;
    }
}

const float GammaFactor = 2.2f;

void Image::gammaCorrection() {
   auto lumBound = [](float Clum) -> float {
        return tanh(3.0f * (Clum-0.75f)) / 3.0f + 0.75f;
    };
    for (int i = 0; i < width * height; ++i) {
        pixelarray[i] = glm::max(pixelarray[i], 0.0f);
        float Clum = dot(pixelarray[i], RGB_Weight);
        if (Clum > 0.75f)
            pixelarray[i] = pixelarray[i] / Clum * lumBound(Clum);
        pixelarray[i] = glm::min(pixelarray[i], 1.0f);
        pixelarray[i] = glm::pow(pixelarray[i], vec3(1.0f / GammaFactor));
    }
}

void Image::postProcessing(int shadeOptions, float exposure) {
    shade(exposure, shadeOptions);
    if (shadeOptions & DoBloom)
        bloom();
    if (shadeOptions & DoFXAA)
        FXAA();
    gammaCorrection();
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
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            int id = y * width + x;
            png_byte *row = &image_data[id * 3];
            for (int k = 0; k < 3; ++k)
                row[k] = static_cast<png_byte>(pixelarray[id][k] * 255);
        }

    for (int y = 0; y < height; ++y)
        png_write_row(png_ptr, &image_data[y * width * 3]);

    png_write_end(png_ptr, nullptr);

    fclose(fp);
    png_destroy_write_struct(&png_ptr, &info_ptr);
}

void Image::load(const char* file_name) {
    FILE *fp = fopen(file_name, "rb");
    if (!fp) {
        std::cerr << "Could not open file for reading" << std::endl;
        return;
    }

    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png_ptr) {
        std::cerr << "Could not create read struct" << std::endl;
        fclose(fp);
        return;
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        std::cerr << "Could not create info struct" << std::endl;
        png_destroy_read_struct(&png_ptr, nullptr, nullptr);
        fclose(fp);
        return;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
        std::cerr << "Error during png creation" << std::endl;
        png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
        fclose(fp);
        return;
    }

    png_init_io(png_ptr, fp);
    png_read_info(png_ptr, info_ptr);

    width = png_get_image_width(png_ptr, info_ptr);
    height = png_get_image_height(png_ptr, info_ptr);
    png_byte color_type = png_get_color_type(png_ptr, info_ptr);
    png_byte bit_depth = png_get_bit_depth(png_ptr, info_ptr);

    std::cout << "width: " << width << " height: " << height << std::endl;
    std::cout << "color_type: " << (int)color_type << " bit_depth: " << (int)bit_depth << std::endl;

    if (bit_depth == 16)
        png_set_strip_16(png_ptr);

    if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png_ptr);

    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
        png_set_expand_gray_1_2_4_to_8(png_ptr);

    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
        png_set_tRNS_to_alpha(png_ptr);

    if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_filler(png_ptr, 0xFF, PNG_FILLER_AFTER);

    if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(png_ptr);

    png_read_update_info(png_ptr, info_ptr);

    pixelarray = new vec3[width * height];

    std::vector<png_byte> image_data(width * height * 4);
    for (int y = 0; y < height; ++y) {
        png_bytep row = &image_data[y * width * 4];
        png_read_row(png_ptr, row, nullptr);
        for (int x = 0; x < width; ++x) {
            int id = y * width + x;
            pixelarray[id] = vec3(row[x * 4] / 255.0f, row[x * 4 + 1] / 255.0f, row[x * 4 + 2] / 255.0f);
        }
    }

    fclose(fp);
    png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
}

void Image::reverseGammaCorrection() {
    for (int i = 0; i < width * height; ++i) {
        pixelarray[i] = glm::pow(pixelarray[i], vec3(GammaFactor));
    }
}

void accumulateInwardRadiance_basic(RadianceData &radiance, vec3 inradiance, float weight) {
    if (!isfinite(inradiance)) {
        std::cerr << "NaN detected!" << std::endl;
        throw std::runtime_error("NaN detected!");
    }
    if (!isfinite(weight)) {
        std::cerr << "NaN detected! (weight)" << std::endl;
        throw std::runtime_error("NaN detected!");
    }
    radiance.radiance += inradiance * weight;
    radiance.Var += dot(inradiance, inradiance) * weight;
}

void accumulateInwardRadiance(const vec3 &baseColor, const LightSample &sample,
                              RadianceData &radiance_d, RadianceData &radiance_s) {

    if (length(sample.light) < eps_zero)
        return;
    vec3 baseColor0 = normalize(baseColor);
    static const vec3 White = normalize(vec3(1.0f));
    float XdotY = dot(baseColor0, White);
    // Consider Pdf = a * white + b * baseColor0
    vec3 perp = normalize(cross(baseColor0, White));
    vec3 bsdfPdf_p = sample.bsdfPdf - perp * dot(perp, sample.bsdfPdf);
    float d1 = dot(bsdfPdf_p, White);
    float d2 = dot(bsdfPdf_p, baseColor0);
    float AplusB = (d1+d2) / (1+XdotY);
    float AminusB = (d1-d2) / (1-XdotY);
    float B = (AplusB - AminusB) / 2.0f;
    vec3 brdfPdf_base = B * baseColor0;
    accumulateInwardRadiance_basic(
            radiance_d, sample.light * B / length(baseColor), sample.weight);
    accumulateInwardRadiance_basic(
            radiance_s, sample.light * (sample.bsdfPdf - brdfPdf_base), sample.weight);
}