#include <thread>
#include <random>
#include <iostream>
#include "render.h"
#include "image.h"
#include "sampling.h"

RenderData::RenderData(int seed) : gen(seed) {
    T_RayAndTexture = C_lightSamples = 0;
}

const float P_RR = 0.6f;
const int maxRayDepth = 16;

vec3 sampleRay(Ray ray, const Model &model, RenderData &renderData, int depth) {

    renderData.T_RayAndTexture -= clock();
    HitInfo info;
    model.rayHit(ray, info);
    renderData.T_RayAndTexture += clock();

    if (info.t == INFINITY)
        return vec3(0.0f);

    if (length(info.emission) > 0)
        return vec3(0.0f); // 不再计算直接光照

    vec3 intersection = ray.origin + ray.direction * info.t;
    BRDF brdf(ray.direction, info.shapeNormal, info.surfaceNormal);

    if (depth == maxRayDepth || renderData.gen() > P_RR) {
        renderData.C_lightSamples++;
        vec3 light = sampleDirectLight(intersection, brdf, model, renderData.gen);
        return (vec3)info.diffuseColor * light / ((depth < maxRayDepth) ? (1.0f - P_RR) : 1.0f);
        // RR 停止时，用光源重要性采样计算直接光照
    }

    float sampleFactor;
    vec3 newDirection = brdf.sample(renderData.gen, sampleFactor);
    if (sampleFactor == 0.0f)
        return vec3(0.0f);
    ray = {intersection, newDirection};
    vec3 inradiance = sampleRay(ray, model, renderData, depth + 1) * sampleFactor / P_RR;
    return (vec3)info.diffuseColor * inradiance;
}

vec3 sampleDirectLightFromFirstIntersection(const GbufferData &Gbuffer, const vec3 &inDir, const vec3 &emission, const Model &model, RenderData &renderData) {

    if (Gbuffer.face == nullptr || length(emission) > 0)
        return vec3(0.0f); // 暂不计入直接光照

    BRDF brdf(inDir, Gbuffer.shapeNormal, Gbuffer.surfaceNormal);
    return sampleDirectLight(Gbuffer.position, brdf, model, renderData.gen);
}

vec3 sampleRayFromFirstIntersection(const GbufferData &Gbuffer, const vec3 &inDir, const vec3 &emission, const Model &model, RenderData &renderData) {

    if (Gbuffer.face == nullptr || length(emission) > 0)
        return vec3(0.0f); // 暂不计入直接光照

    BRDF brdf(inDir, Gbuffer.shapeNormal, Gbuffer.surfaceNormal);

    float sampleFactor;
    vec3 newDirection = brdf.sample(renderData.gen, sampleFactor);
    if (sampleFactor == 0.0f)
        return vec3(0.0f);
    Ray ray = {Gbuffer.position, newDirection};
    return sampleRay(ray, model, renderData, 1) * sampleFactor;
}

void rayCasting(Model &model, const RenderArgs &args, Image &image) {
    const int
            width = args.width,
            height = args.height;
    vec3
            direction = args.direction,
            right = args.right,
            up = args.up,
            position = args.position;
    const float
            accuracy = args.accuracy;
    for (int x = 0; x < width; x++)
        for (int y = 0; y < height; y++) {
            float
                    rayX = x - width / 2.0f,
                    rayY = y - height / 2.0f;
            vec3 aim = normalize(direction + accuracy * (rayX * right + rayY * up));
            Ray ray = {position, aim};
            GbufferData &Gbuffer = image.Gbuffer[y * width + x];
            float depth = INFINITY;
            model.rayCast(ray, depth, Gbuffer.face);

            if (Gbuffer.face == nullptr)
                continue;

            vec3 intersection = position + aim * depth;
            Gbuffer.position = intersection;
            HitInfo hitInfo;
            getHitInfo(*Gbuffer.face, intersection, aim, hitInfo);
            Gbuffer.shapeNormal = hitInfo.shapeNormal;
            Gbuffer.surfaceNormal = hitInfo.surfaceNormal;
        }
}

void normalize(RadianceData &radiance, int spp) {
    radiance.radiance /= (float) spp;
    radiance.Var /= (float) spp;
    radiance.Var -= dot(radiance.radiance, radiance.radiance);
    if (radiance.Var < 0.0f)
        radiance.Var = 0.0f;
}

void renderPixel(const Model &model, const RenderArgs &args,
                 RenderData &renderData, Image &image, int id) {
    GbufferData &Gbuffer = image.Gbuffer[id];
    if (Gbuffer.face == nullptr)
        return;
    const int
            spp_direct = args.spp * args.P_Direct,
            spp_indirect = args.spp - spp_direct;
    const float
            exposure = args.exposure;
    vec3 aim = normalize(Gbuffer.position - args.position);
    RadianceData &radiance_d = image.radiance_d[id];
    RadianceData &radiance_i = image.radiance_i[id];
    HitInfo hitInfo;
    getHitInfo(*Gbuffer.face, Gbuffer.position, aim, hitInfo);
    for (int T = 0; T < spp_direct; T++) {
        vec3 inradiance_d = sampleDirectLightFromFirstIntersection(Gbuffer, aim, hitInfo.emission, model, renderData);
        inradiance_d *= exposure;
        radiance_d.radiance += inradiance_d;
        radiance_d.Var += dot(inradiance_d, inradiance_d);
    }
    normalize(radiance_d, spp_direct);
    for (int T = 0; T < spp_indirect; T++) {
        vec3 inradiance_i = sampleRayFromFirstIntersection(Gbuffer, aim, hitInfo.emission, model, renderData);
        inradiance_i *= exposure;
        radiance_i.radiance += inradiance_i;
        radiance_i.Var += dot(inradiance_i, inradiance_i);
    }
    normalize(radiance_i, spp_indirect);
}

void render(const Model &model, const RenderArgs &args,
            RenderData &renderData, Image &image, int xL, int xR) {
    static const int
            C_calc = 256,
            C_report = 4096;
    const int
            width = args.width,
            height = args.height;
    int cnt = 0;
    for (int y = 0; y < height; y++)
        for (int x = xL; x < xR; x++) {
            renderPixel(model, args, renderData, image, y * args.width + x);
            cnt++;
            if (cnt == C_calc) {
                cnt = 0;
                (*renderData.renderedPixels) += C_calc;
                if ((*renderData.renderedPixels) % C_report == 0)
                    std::cout << "Rendered " << 1e2f *(*renderData.renderedPixels)/(width * height) << "%" << std::endl;
            }
        }
}

void render_multiThread(Model &model, const RenderArgs &args) {
    const int
            startTime = clock(),
            width = args.width,
            height = args.height,
            threads = args.threads;

    std::cout << "Rendering started with " << threads << " threads." << std::endl;

    Image image(width, height);
    rayCasting(model, args, image);

    std::cout << "Ray casting completed. (" << clock()-startTime << " ms)"<< std::endl;

    image.shade(args.position, Image::DiffuseColor);
    image.save((args.savePath+"(DiffuseColor).png").c_str());

    std::atomic<int> renderedPixels(0);
    std::vector<RenderData> datas;
    datas.reserve(threads);
    for (int i = 0; i < threads; i++) {
        datas.emplace_back(i);
        datas.back().renderedPixels = &renderedPixels;
    }

    std::vector<std::thread> threadPool;
    threadPool.reserve(threads);
    for (int i = 0; i < threads; i++)
        threadPool.emplace_back([&](int index) {
            render(model, args, datas[index], image,
                     index * width / threads, (index + 1) * width / threads);
        }, i);
    for (auto &thread : threadPool)
        thread.join();

    int T_ray = 0, C_lightSamples = 0;
    for (const auto &data : datas) {
        T_ray += data.T_RayAndTexture;
        C_lightSamples += data.C_lightSamples;
    }

    std::cout << "Rendering completed in " << clock() - startTime << " ms." << std::endl;
    std::cout << "Ray intersection & Texture query time: " << T_ray << " ms*thread." << std::endl;
    std::cout << "Direct light samples: " << C_lightSamples << std::endl;

    image.shade(args.position, Image::DirectLight);
    image.save((args.savePath+"(DirectLight).png").c_str());
    image.shade(args.position, Image::IndirectLight);
    image.save((args.savePath+"(IndirectLight).png").c_str());
    image.shade(args.position, Image::DirectLight | Image::IndirectLight | Image::Emission | Image::DiffuseColor);
    image.save((args.savePath+"(raw).png").c_str());
    image.filter();
    image.shade(args.position, Image::DirectLight);
    image.save((args.savePath+"(DirectLight.filtered).png").c_str());
    image.shade(args.position, Image::IndirectLight);
    image.save((args.savePath+"(IndirectLight.filtered).png").c_str());
    image.shade(args.position, Image::DirectLight | Image::IndirectLight | Image::Emission | Image::DiffuseColor);
    image.save((args.savePath+"(filtered).png").c_str());

    std::cout << "Post processing finished. Total: " << clock() - startTime << " ms." << std::endl;
}