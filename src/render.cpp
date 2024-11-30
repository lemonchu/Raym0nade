#include <thread>
#include <random>
#include <iostream>
#include "render.h"
#include "image.h"
#include "sampling.h"

const float P_RR = 0.75f;
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

    if (depth == maxRayDepth || std::uniform_real_distribution<float>(0.0f, 1.0f)(renderData.gen) > P_RR) {
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
    const unsigned int
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
    for (unsigned int x = 0; x < width; x++)
        for (unsigned int y = 0; y < height; y++) {
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
                getHitInfo(*Gbuffer.face, intersection, hitInfo);
                Gbuffer.position = intersection;
                Gbuffer.diffuseColor = (vec3)hitInfo.diffuseColor;
                Gbuffer.shapeNormal = hitInfo.shapeNormal;
                Gbuffer.surfaceNormal = hitInfo.surfaceNormal;
                Gbuffer.emission = hitInfo.emission * exposure;
            }
        }
}

const float P_Direct = 0.25f;

void render(const Model &model, const RenderArgs &args, RenderData &renderData, Image &image, int xL, int xR) {
    const unsigned int
            width = args.width,
            height = args.height,
            spp_direct = args.spp * P_Direct,
            spp_indirect = args.spp - spp_direct,
            spp = args.spp;
    const float
            accuracy = args.accuracy,
            exposure = args.exposure;
    vec3
            direction = args.direction,
            right = args.right,
            up = args.up,
            position = args.position;

    for (int T = 0; T < spp_direct; T++) {
        for (unsigned int x = xL; x < xR; x++)
            for (unsigned int y = 0; y < height; y++) {
                float
                        rayX = x - width / 2.0f,
                        rayY = y - height / 2.0f;
                vec3 aim = normalize(direction + accuracy * (rayX * right + rayY * up));

                RadianceData &radiance = image.radiance_d[y * width + x];
                GbufferData &Gbuffer = image.Gbuffer[y * width + x];
                vec3 intersection = position + aim * Gbuffer.depth;
                vec3 inradiance = sampleDirectLightFromFirstIntersection(Gbuffer, aim, model, renderData);
                inradiance *= exposure;
                radiance.radiance += inradiance;
                radiance.Var += dot(inradiance, inradiance);
                radiance.sampleCount++;
            }
        std::cout << "Thread " << std::this_thread::get_id() << " progress: " << T + 1 << " / " << spp << std::endl;
    }

    for (int T = 0; T < spp_indirect; T++) {
        for (unsigned int x = xL; x < xR; x++)
            for (unsigned int y = 0; y < height; y++) {
                float
                        rayX = x - width / 2.0f,
                        rayY = y - height / 2.0f;
                vec3 aim = normalize(direction + accuracy * (rayX * right + rayY * up));

                RadianceData &pixel = image.radiance_i[y * width + x];
                GbufferData &Gbuffer = image.Gbuffer[y * width + x];
                vec3 intersection = position + aim * Gbuffer.depth;
                vec3 inradiance = sampleRayFromFirstIntersection(Gbuffer, aim, intersection, model, renderData);
                inradiance *= exposure;
                pixel.radiance += inradiance;
                pixel.Var += dot(inradiance, inradiance);
                pixel.sampleCount++;
            }
        std::cout << "Thread " << std::this_thread::get_id() << " progress: " << T + spp_direct + 1 << " / " << spp << std::endl;
    }
}

RenderData::RenderData(int seed) : gen(seed) {
    T_RayAndTexture = 0;
}

void render_multiThread(Model &model, const RenderArgs &args) {
    const unsigned int
            startTime = clock(),
            width = args.width,
            height = args.height,
            threads = args.threads;

    std::cout << "Rendering started with " << threads << " threads." << std::endl;

    Image image(width, height);
    rayCasting(model, args, image);

    std::cout << "Ray casting completed." << std::endl;

    std::vector<RenderData> datas;
    datas.reserve(threads);
    for (unsigned int i = 0; i < threads; ++i)
        datas.emplace_back(0);
    auto renderTask = [&](int threadIndex) {
        render(model, args, datas[threadIndex], image,
               width/threads * threadIndex, width/threads * (threadIndex+1));
    };

    std::vector<std::thread> threadPool;
    for (unsigned int i = 0; i < threads; ++i)
        threadPool.emplace_back(renderTask, i);
    for (auto &thread : threadPool)
        thread.join();

    std::cout << "Rendering completed in " << clock() - startTime << " ms." << std::endl;

    int T_ray = 0;
    for (auto &data : datas)
        T_ray += data.T_RayAndTexture;

    std::cout << "Ray intersection & Texture query time: " << T_ray << " ms*thread." << std::endl;

    image.normalizeRadiance();
    image.filter();
    image.shade();
    image.bloom();
    image.gammaCorrection();
    image.save(args.savePath.c_str());

    std::cout << "Post processing finished. Total: " << clock() - startTime << " ms." << std::endl;
}