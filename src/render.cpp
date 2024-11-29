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

const float P_Direct = 0.25f;

vec3 sampleRayFromFirstIntersection(const Face* face, const vec3 &inDir, const vec3 &intersection, const Model &model, RenderData &renderData) {

    if (face == nullptr)
        return vec3(0.0f);

    HitInfo info;
    getHitInfo(*face, intersection, info);

    if (length(info.emission) > 0)
        return vec3(0.0f); // 暂不计入直接光照

    BRDF brdf(inDir, info.shapeNormal, info.surfaceNormal);
    //if (std::uniform_real_distribution<float>(0.0f, 1.0f)(renderData.gen) < P_Direct)
        return sampleDirectLight(intersection, brdf, model, renderData.gen);// / P_Direct;

    /*float sampleFactor;
    vec3 newDirection = brdf.sample(renderData.gen, sampleFactor);
    if (sampleFactor == 0.0f)
        return vec3(0.0f);
    Ray ray = {intersection, newDirection};
    return sampleRay(ray, model, renderData, 1) * sampleFactor / (1.0f - P_Direct);*/
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
            PixelData &pixel = image.buffer[y * width + x];
            model.rayCast(ray, pixel.depth, pixel.face);

            if (pixel.face == nullptr)
                pixel.diffuseColor = pixel.shapeNormal = vec3(0.0f);
            else {
                vec3 intersection = position + aim * pixel.depth;
                HitInfo hitInfo;
                getHitInfo(*pixel.face, intersection, hitInfo);
                pixel.position = intersection;
                pixel.diffuseColor = (vec3)hitInfo.diffuseColor;
                pixel.shapeNormal = hitInfo.shapeNormal;
                pixel.surfaceNormal = hitInfo.surfaceNormal;
                pixel.emission = hitInfo.emission * exposure;
            }
        }
}

void render(const Model &model, const RenderArgs &args, RenderData &renderData, Image &image, int xL, int xR) {
    const unsigned int
            width = args.width,
            height = args.height,
            spp = args.spp;
    const float
            accuracy = args.accuracy,
            exposure = args.exposure;
    vec3
            direction = args.direction,
            right = args.right,
            up = args.up,
            position = args.position;
    for (int T = 0; T < spp; T++) {
        for (unsigned int x = xL; x < xR; x++)
            for (unsigned int y = 0; y < height; y++) {
                float
                        rayX = x - width / 2.0f,
                        rayY = y - height / 2.0f;
                vec3 aim = normalize(direction + accuracy * (rayX * right + rayY * up));

                PixelData &pixel = image.buffer[y * width + x];
                vec3 intersection = position + aim * pixel.depth;
                vec3 inradiance = sampleRayFromFirstIntersection(pixel.face, aim, intersection, model, renderData) * exposure;
                pixel.inradiance += inradiance;
                pixel.Epower2 += dot(inradiance, inradiance);
                pixel.sampleCount++;
            }
        std::cout << "Thread " << std::this_thread::get_id() << " progress: " << T + 1 << " / " << spp << std::endl;
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

    image.filterInradiance();

    image.save(args.savePath.c_str());

    std::cout << "Rendering completed in " << clock() - startTime << " ms." << std::endl;

    int T_ray = 0;
    for (auto &data : datas)
        T_ray += data.T_RayAndTexture;
    std::cout << "Ray intersection & Texture query time: " << T_ray << " ms*thread." << std::endl;
}