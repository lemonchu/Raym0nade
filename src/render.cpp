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

vec3 sampleDirectLightFromFirstIntersection(const GbufferData &Gbuffer, const vec3 &inDir, const Model &model, RenderData &renderData) {

    if (Gbuffer.face == nullptr || length(Gbuffer.emission) > 0)
        return vec3(0.0f); // 暂不计入直接光照

    BRDF brdf(inDir, Gbuffer.shapeNormal, Gbuffer.surfaceNormal);
    return sampleDirectLight(Gbuffer.position, brdf, model, renderData.gen);
}

vec3 sampleRayFromFirstIntersection(const GbufferData &Gbuffer, const vec3 &inDir, const vec3 &intersection, const Model &model, RenderData &renderData) {

    if (Gbuffer.face == nullptr || length(Gbuffer.emission) > 0)
        return vec3(0.0f); // 暂不计入直接光照

    BRDF brdf(inDir, Gbuffer.shapeNormal, Gbuffer.surfaceNormal);

    float sampleFactor;
    vec3 newDirection = brdf.sample(renderData.gen, sampleFactor);
    if (sampleFactor == 0.0f)
        return vec3(0.0f);
    Ray ray = {intersection, newDirection};
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
            accuracy = args.accuracy,
            exposure = args.exposure;
    for (int x = 0; x < width; x++)
        for (int y = 0; y < height; y++) {
            float
                    rayX = x - width / 2.0f,
                    rayY = y - height / 2.0f;
            vec3 aim = normalize(direction + accuracy * (rayX * right + rayY * up));
            Ray ray = {position, aim};
            GbufferData &Gbuffer = image.Gbuffer[y * width + x];
            model.rayCast(ray, Gbuffer.depth, Gbuffer.face);

            if (Gbuffer.face == nullptr)
                Gbuffer.diffuseColor = Gbuffer.shapeNormal = vec3(0.0f);
            else {
                vec3 intersection = position + aim * Gbuffer.depth;
                HitInfo hitInfo;
                getHitInfo(*Gbuffer.face, intersection, aim, hitInfo);
                Gbuffer.position = intersection;
                Gbuffer.diffuseColor = (vec3)hitInfo.diffuseColor;
                Gbuffer.shapeNormal = hitInfo.shapeNormal;
                Gbuffer.surfaceNormal = hitInfo.surfaceNormal;
                Gbuffer.emission = hitInfo.emission * exposure;
            }
        }
}

void render_d(const Model &model, const RenderArgs &args,
            RenderData &renderData, Image &image, int xL, int xR) {
    const int
            width = args.width,
            height = args.height;
    const float
            exposure = args.exposure,
            accuracy = args.accuracy;
    const vec3
            position = args.position,
            direction = args.direction,
            right = args.right,
            up = args.up;
    for (int x = xL; x < xR; x++)
        for (int y = 0; y < args.height; y++) {
            float
                    rayX = x - width / 2.0f,
                    rayY = y - height / 2.0f;
            vec3 aim = normalize(direction + accuracy * (rayX * right + rayY * up));
            RadianceData &radiance_d = image.radiance_d[y * width + x];
            GbufferData &Gbuffer = image.Gbuffer[y * width + x];
            vec3 inradiance = sampleDirectLightFromFirstIntersection(Gbuffer, aim, model, renderData);
            inradiance *= exposure;
            radiance_d.radiance += inradiance;
            radiance_d.Spower2 += dot(inradiance, inradiance);
            radiance_d.sampleCount++;
        }
}

void render_i(const Model &model, const RenderArgs &args,
            RenderData &renderData, Image &image, int xL, int xR) {
    const int
            width = args.width,
            height = args.height;
    const float
            exposure = args.exposure,
            accuracy = args.accuracy;
    const vec3
            position = args.position,
            direction = args.direction,
            right = args.right,
            up = args.up;
    for (int x = xL; x < xR; x++)
        for (int y = 0; y < args.height; y++) {
            float
                    rayX = x - width / 2.0f,
                    rayY = y - height / 2.0f;
            vec3 aim = normalize(direction + accuracy * (rayX * right + rayY * up));
            RadianceData &radiance_i = image.radiance_i[y * width + x];
            GbufferData &Gbuffer = image.Gbuffer[y * width + x];
            vec3 inradiance = sampleRayFromFirstIntersection(Gbuffer, aim, Gbuffer.position, model, renderData);
            inradiance *= exposure;
            radiance_i.radiance += inradiance;
            radiance_i.Spower2 += dot(inradiance, inradiance);
            radiance_i.sampleCount++;
        }
}

void render_multiThread(Model &model, const RenderArgs &args) {
    const int
            startTime = clock(),
            width = args.width,
            height = args.height,
            threads = args.threads,
            spp_direct = args.spp * args.P_Direct,
            spp_indirect = args.spp - spp_direct;

    std::cout << "Rendering started with " << threads << " threads." << std::endl;
    int sav = clock();

    Image image(width, height);
    rayCasting(model, args, image);

    std::cout << "Ray casting completed. (" << clock()-sav << " ms)"<< std::endl;
    sav = clock();

    image.shade(Image::DiffuseColor);
    image.save((args.savePath+"(DiffuseColor).png").c_str());

    std::vector<RenderData> datas;
    datas.reserve(threads);
    for (int i = 0; i < threads; i++)
        datas.emplace_back(i);

    for (int T = 0; T < spp_direct; T++) {
        std::vector<std::thread> threadPool;
        threadPool.resize(threads);
        for (int i = 0; i < threads; i++)
            threadPool.emplace_back([&](int index) {
                render_d(model, args, datas[index], image,
                         index * width / threads, (index + 1) * width / threads);
            }, i);
        for (auto &thread : threadPool)
            thread.join();

        std::cout << "Direct light rendering: " << 1e2f*(T+1)/ spp_direct << "% (" << clock() - sav << " ms)" << std::endl;
    }

    sav = clock();

    for (int T = 0; T < spp_indirect; T++) {
        std::vector<std::thread> threadPool;
        threadPool.resize(threads);
        for (int i = 0; i < threads; i++)
            threadPool.emplace_back([&](int index) {
                render_i(model, args, datas[index], image,
                         index * width / threads, (index + 1) * width / threads);
            }, i);
        for (auto &thread : threadPool)
            thread.join();

        std::cout << "Indirect light rendering: " << 1e2f*(T+1)/spp_indirect << "% (" << clock() - sav << " ms)" << std::endl;
    }

    int T_ray = 0, C_lightSamples = 0;
    for (auto &data : datas) {
        T_ray += data.T_RayAndTexture;
        C_lightSamples += data.C_lightSamples;
    }

    std::cout << "Rendering completed in " << clock() - startTime << " ms." << std::endl;
    std::cout << "Ray intersection & Texture query time: " << T_ray << " ms*thread." << std::endl;
    std::cout << "Direct light samples: " << C_lightSamples << std::endl;

    image.normalizeRadiance();
    image.shade(Image::DirectLight | Image::ShowVar);
    image.save((args.savePath+"(DirectLight.Var).png").c_str());
    image.shade(Image::IndirectLight | Image::ShowVar);
    image.save((args.savePath+"(IndirectLight.Var).png").c_str());
    image.shade(Image::DirectLight);
    image.save((args.savePath+"(DirectLight).png").c_str());
    image.shade(Image::IndirectLight);
    image.save((args.savePath+"(IndirectLight).png").c_str());
    image.shade(Image::DirectLight | Image::IndirectLight | Image::Emission | Image::DiffuseColor);
    image.save((args.savePath+"(raw).png").c_str());
    image.filter();
    image.shade(Image::DirectLight);
    image.save((args.savePath+"(DirectLight.filtered).png").c_str());
    image.shade(Image::IndirectLight);
    image.save((args.savePath+"(IndirectLight.filtered).png").c_str());
    image.shade(Image::DirectLight | Image::IndirectLight | Image::Emission | Image::DiffuseColor);
    image.save((args.savePath+"(filtered).png").c_str());

    std::cout << "Post processing finished. Total: " << clock() - startTime << " ms." << std::endl;
}