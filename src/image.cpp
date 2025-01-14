#include <vector>
#include <png.h>
#include <iostream>
#include <thread>
#include "image.h"

RadianceData::RadianceData() : radiance(vec3(0.0f)), Var(0.0f) {}

Photo::Photo(int width, int height) : width(width), height(height) {
    Gbuffer = new HitInfo[width * height];
    radiance_Dd = new RadianceData[width * height];
    radiance_Ds = new RadianceData[width * height];
    radiance_Id = new RadianceData[width * height];
    radiance_Is = new RadianceData[width * height];
    pixelarray = new vec3[width * height];
}

Photo::Photo(const char *file_name) : radiance_Dd(nullptr), radiance_Ds(nullptr), radiance_Id(nullptr), radiance_Is(nullptr), pixelarray(nullptr), Gbuffer(nullptr) {
    load(file_name);
}

Photo::~Photo() {
    delete[] Gbuffer;
    delete[] radiance_Dd;
    delete[] radiance_Ds;
    delete[] radiance_Id;
    delete[] radiance_Is;
}

void clamp(RadianceData *radiance, int width, int height) {
    const int filterRadius = 3;
    const float h[7] = {0.03125, 0.109375, 0.21875, 0.28125, 0.21875, 0.109375, 0.03125};
    std::vector<float> Clum(width * height, 0.0f);
    std::vector<float> blurredClum1(width * height, 0.0f);
    std::vector<float> blurredClum2(width * height, 0.0f);

    // Calculate luminance for each pixel
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            Clum[y * width + x] = dot(radiance[y * width + x].radiance, RGB_Weight);

    // Apply horizontal 7x1 Gaussian blur to the luminance values
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            float sum = 0.0f;
            for (int dx = -filterRadius; dx <= filterRadius; ++dx) {
                int nx = x + dx;
                if (nx >= 0 && nx < width)
                    sum += Clum[y * width + nx] * h[dx + filterRadius];
            }
            blurredClum1[y * width + x] = sum;
        }

    // Apply vertical 1x7 Gaussian blur to the luminance values
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            float sum = 0.0f;
            for (int dy = -filterRadius; dy <= filterRadius; ++dy) {
                int ny = y + dy;
                if (ny >= 0 && ny < height)
                    sum += blurredClum1[ny * width + x] * h[dy + filterRadius];
            }
            blurredClum2[y * width + x] = sum;
        }

    static const float clampThreshold = 36.0f;
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            float &lum = Clum[y * width + x];
            float blurredLum =
                    blurredClum2[y * width + x] - h[filterRadius] * h[filterRadius] * lum;
            if (lum > clampThreshold * blurredLum)
                radiance[y * width + x].radiance *= blurredLum / (lum + eps_zero) / (1.0f - h[filterRadius] * h[filterRadius]);
        }
}

void Photo::spatialClamp() {
    clamp(radiance_Dd, width, height);
    clamp(radiance_Ds, width, height);
    clamp(radiance_Id, width, height);
    clamp(radiance_Is, width, height);
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

float getWeight(const HitInfo &Gp, const HitInfo &Gq, const RadianceData &Lp, const RadianceData &Lq, bool specular) {
    static constexpr float
            sigma_z = 1.0f,
            sigma_n = 1024.0f,
            sigma_l = 1.0f,
            sigma_m = 1.0f,
            eps = 1e-2f,
            eps_r = 4e-2f;

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

    if (specular)
        weight *= std::max(Gp.roughness, eps_r);

    if (!isfinite(weight))
        return 0.0f;

    return weight;
}

void filterRadiance(RadianceData *radiance, const HitInfo *Gbuffer, int width, int height, bool specular, int step) {

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
                        weight *= getWeight(Gp, Gq, Lp, Lq, specular);
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

void filterRadiance(RadianceData *radiance, const HitInfo *Gbuffer, int width, int height, bool specular) {
    filterVar(radiance, width, height);
    filterRadiance(radiance, Gbuffer, width, height, specular, 1);
    filterRadiance(radiance, Gbuffer, width, height, specular,2);
    filterRadiance(radiance, Gbuffer, width, height, specular,4);
    filterRadiance(radiance, Gbuffer, width, height, specular,8);
    filterRadiance(radiance, Gbuffer, width, height, specular,16);
}

void Photo::filter() {
    std::thread t1([this]() { filterRadiance(radiance_Dd, Gbuffer, width, height, false); });
    std::thread t2([this]() { filterRadiance(radiance_Ds, Gbuffer, width, height, true); });
    std::thread t3([this]() { filterRadiance(radiance_Id, Gbuffer, width, height, false); });
    std::thread t4([this]() { filterRadiance(radiance_Is, Gbuffer, width, height, true); });

    t1.join();
    t2.join();
    t3.join();
    t4.join();
}

void Photo::shade(int options) {
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

void Photo::bloom() {
    static const int filterRadius = 2;
    static const float h[5] = {0.0625, 0.25, 0.375, 0.25, 0.0625};
    std::vector<vec3> bloomColor(width * height);

    static const float omega_bloom = 0.65f; // Bloom 强度
    for (int id = 0; id < width * height; ++id) {
        float Clum = dot(pixelarray[id], RGB_Weight);
        if (Clum < 1.0f)
            continue;
        bloomColor[id] = pixelarray[id] / pow(Clum, omega_bloom);
        for (int k = 0; k < 3; ++k)
            bloomColor[id][k] = std::max(bloomColor[id][k] - 1.0f, 0.0f);
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
            pixelarray[id] += bloomColor[id] / 6.0f;
    }
}

void Photo::depthFeildBlur() {
    struct PixelPos {
        int x, y;
        float depth;
        PixelPos(int x, int y, float depth) : x(x), y(y), depth(depth) {}
    };
    std::vector<PixelPos> pixels;
    std::vector<vec3> blurredPixelArray(width * height, vec3(0.0f));
    std::vector<float> gained(width * height, 0.0f);

    pixels.reserve(width * height);
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            float depth = length(Gbuffer[y * width + x].position - cameraPosition);
            pixels.emplace_back(x, y, depth);
        }

    auto cmp = [](const PixelPos &a, const PixelPos &b) {
        return a.depth < b.depth;
    };
    std::stable_sort(pixels.begin(), pixels.end(), cmp);

    constexpr float gainThreshold = 0.99f;
    for (int i = 0; i < width*height; i++) {
        const PixelPos &pixel = pixels[i];
        int x = pixel.x, y = pixel.y;

        constexpr float MaxCoC = 96.0f;
        float CoC0 = std::min(CoC * abs(1.0f - focus / pixel.depth) + eps_zero, MaxCoC);

        int radius = int(CoC0);
        vec3 color = pixelarray[y * width + x];
        float totalWeight = 0.0f;

        for (int dy = -radius; dy <= radius; ++dy) {
            int xLen = int(sqrt(CoC0 * CoC0 - dy * dy));
            int xl, xr;
            for (xl = -xLen; xl <= 0; xl++) {
                float weight = std::min(CoC0 - sqrtf(xl * xl + dy * dy), 1.0f);
                if (weight == 1.0f)
                    break;
                totalWeight += weight;
            }
            for (xr = xLen; xr > 0; xr--) {
                float weight = std::min(CoC0 - sqrtf(xr * xr + dy * dy), 1.0f);
                if (weight == 1.0f)
                    break;
                totalWeight += weight;
            }
            totalWeight += xr - xl + 1;
        }

        for (int dy = -radius; dy <= radius; ++dy) {
            int xLen = int(sqrt(CoC0 * CoC0 - dy * dy));
            for (int dx = -xLen; dx <= xLen; ++dx) {
                int nx = x + dx, ny = y + dy;
                if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                    continue;
                float weight = std::min(CoC0 - sqrtf(dx * dx + dy * dy), 1.0f) / totalWeight;
                if (weight < eps_zero)
                    continue;
                int id = ny * width + nx;
                if (gained[id] + weight > gainThreshold)
                    weight = gainThreshold - gained[id];
                if (weight < eps_zero)
                    continue;
                blurredPixelArray[id] += color * weight;
                gained[id] += weight;
            }
        }
    }

    for (int i = 0; i < width * height; ++i)
        pixelarray[i] = blurredPixelArray[i] / gained[i];
}

const float EDGE_THRESHOLD_MIN = 0.0312f;
const float EDGE_THRESHOLD_MAX = 0.125f;
const float SUBPIXEL_QUALITY = 0.75f;
const int QUALITY = 12;

void Photo::FXAA() {

    const vec3 *data = pixelarray;
    vec3 *output = new vec3[width * height];

    auto getLuminance = [](const vec3& rgb) -> float {
        return dot(rgb, RGB_Weight);
    };

    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {
            float lumaM = getLuminance(data[y * width + x]);
            float lumaN = y > 0 ? getLuminance(data[(y - 1) * width + x]) : lumaM;
            float lumaS = y < height - 1 ? getLuminance(data[(y + 1) * width + x]) : lumaM;
            float lumaE = x < width - 1 ? getLuminance(data[y * width + x + 1]) : lumaM;
            float lumaW = x > 0 ? getLuminance(data[y * width + x - 1]) : lumaM;

            float rangeMin = std::min({lumaN, lumaS, lumaE, lumaW});
            float rangeMax = std::max({lumaN, lumaS, lumaE, lumaW});
            float range = rangeMax - rangeMin;

            if (range < std::max(EDGE_THRESHOLD_MIN, rangeMax * EDGE_THRESHOLD_MAX)) {
                output[y * width + x] = data[y * width + x];
                continue;
            }

            float lumaNW = (y > 0 && x > 0) ? getLuminance(data[(y - 1) * width + x - 1]) : lumaM;
            float lumaNE = (y > 0 && x < width - 1) ? getLuminance(data[(y - 1) * width + x + 1]) : lumaM;
            float lumaSW = (y < height - 1 && x > 0) ? getLuminance(data[(y + 1) * width + x - 1]) : lumaM;
            float lumaSE = (y < height - 1 && x < width - 1) ? getLuminance(data[(y + 1) * width + x + 1]) : lumaM;

            float edgeHorz = std::abs((lumaNW + lumaW + lumaSW) - (lumaNE + lumaE + lumaSE)) * (1.0f / 3.0f);
            float edgeVert = std::abs((lumaNW + lumaN + lumaNE) - (lumaSW + lumaS + lumaSE)) * (1.0f / 3.0f);

            bool isHorizontal = edgeHorz >= edgeVert;

            float stepLength = isHorizontal ? 1.0f / float(width) : 1.0f / float(height);
            float gradientStep = std::clamp(
                    (isHorizontal ? edgeHorz : edgeVert) / range,
                    -2.0f,
                    2.0f
            );

            vec2 uv(float(x) / float(width), float(y) / float(height));
            vec3 finalColor = data[y * width + x];
            float bestDelta = 0.0f;

            for (int i = 0; i < QUALITY; i++) {
                vec2 offset(
                        isHorizontal ? 0.0f : gradientStep * stepLength * float(i + 1),
                        isHorizontal ? gradientStep * stepLength * float(i + 1) : 0.0f
                );

                vec2 sampleUv = uv + offset;
                if (sampleUv.x < 0.0f || sampleUv.x > 1.0f || sampleUv.y < 0.0f || sampleUv.y > 1.0f) {
                    continue;
                }

                int sx = int(sampleUv.x * float(width));
                int sy = int(sampleUv.y * float(height));
                sx = std::clamp(sx, 0, width - 1);
                sy = std::clamp(sy, 0, height - 1);

                vec3 sampleColor = data[sy * width + sx];
                float sampleLuma = getLuminance(sampleColor);
                float delta = std::abs(sampleLuma - lumaM);

                if (delta > bestDelta) {
                    bestDelta = delta;
                    finalColor = sampleColor;
                }
            }

            float subPixelOffset = std::min(
                    (std::abs(lumaN + lumaS - 2.0f * lumaM) * 2.0f +
                     std::abs(lumaE + lumaW - 2.0f * lumaM)) * 0.25f,
                    1.0f
            );

            output[y * width + x] = mix(
                    data[y * width + x],
                    finalColor,
                    subPixelOffset * SUBPIXEL_QUALITY
            );
        }
    }

    std::memcpy(pixelarray, output, width * height * sizeof(vec3));
    delete[] output;
}

const float GammaFactor = 2.2f;

void Photo::gammaCorrection() {
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

void Photo::postProcessing(int shadeOptions) {
    shade(shadeOptions);
    if (shadeOptions & DoDepthFieldBlur)
        depthFeildBlur();
    if (shadeOptions & DoBloom)
        bloom();
    gammaCorrection();
    if (shadeOptions & DoFXAA)
        FXAA();
}

void Photo::save(const char *file_name) {

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

void Photo::load(const char* file_name) {
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

void Photo::reverseGammaCorrection() {
    for (int i = 0; i < width * height; ++i)
        pixelarray[i] = glm::pow(pixelarray[i], vec3(GammaFactor));
}

void accumulateInwardRadiance_basic(RadianceData &radiance, vec3 inradiance, float weight) {
    if (!isfinite(inradiance)) {
        std::cerr << "NaN detected!" << std::endl;
        //throw std::runtime_error("NaN detected!");
        return;
    }
    if (!isfinite(weight)) {
        std::cerr << "NaN detected! (weight)" << std::endl;
        //throw std::runtime_error("NaN detected!");
        return;
    }
    radiance.radiance += inradiance * weight;
    radiance.Var += dot(inradiance, inradiance) * weight;
}

void accumulateInwardRadiance(const vec3 &baseColor, const LightSample &sample,
                              RadianceData &radiance_d, RadianceData &radiance_s) {

    if (length(sample.light) < eps_zero)
        return;
    vec3 baseColor0 = normalize(baseColor);
    if (length(baseColor) < eps_zero) {
        accumulateInwardRadiance_basic(radiance_s, sample.light * sample.bsdfPdf, sample.weight);
        return;
    }
    static const vec3 White = normalize(vec3(1.0f));
    float XdotY = dot(baseColor0, White);
    if (XdotY > 0.99f) {
        accumulateInwardRadiance_basic(radiance_d, sample.light * sample.bsdfPdf / baseColor, sample.weight);
        return;
    }
    // Consider Pdf = a * white + b * baseColor0
    vec3 perp = normalize(cross(baseColor0, White));
    vec3 bsdfPdf_p = sample.bsdfPdf - perp * dot(perp, sample.bsdfPdf);
    float d1 = dot(bsdfPdf_p, White);
    float d2 = dot(bsdfPdf_p, baseColor0);
    float AplusB = (d1+d2) / (1+XdotY);
    float AminusB = (d1-d2) / (1-XdotY);
    float B = (AplusB - AminusB) / 2.0f;
    vec3 brdfPdf_base = B * baseColor0;
    accumulateInwardRadiance_basic
        (radiance_d, sample.light * B / length(baseColor), sample.weight);
    accumulateInwardRadiance_basic
        (radiance_s, sample.light * (sample.bsdfPdf - brdfPdf_base), sample.weight);
}