#include <thread>
#include <random>
#include <iostream>
#include "render.h"
#include "image.h"
#include "sampling.h"

RenderData::RenderData() {
    T_RayAndTexture = 0;
}

Generator RenderData::generator(int dimension) {
    return {&sobel, dimension};
}

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

    if (depth == maxRayDepth || renderData.sobel() > P_RR) {
        vec3 light = sampleDirectLight(
                intersection, brdf, model,
                renderData.generator(depth*2),
                renderData.generator(depth*2+1),
                renderData.generator()
            );
        return (vec3)info.diffuseColor * light / ((depth < maxRayDepth) ? (1.0f - P_RR) : 1.0f);
        // RR 停止时，用光源重要性采样计算直接光照
    }

    float sampleFactor;
    vec3 newDirection = brdf.sample(
            renderData.generator(depth * 2),
            renderData.generator(depth * 2 + 1),
            sampleFactor
    );
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
    return sampleDirectLight(Gbuffer.position, brdf, model,
             renderData.generator(0),
             renderData.generator(1),
             renderData.generator()
    );
}

vec3 sampleRayFromFirstIntersection(const GbufferData &Gbuffer, const vec3 &inDir, const vec3 &intersection, const Model &model, RenderData &renderData) {

    if (Gbuffer.face == nullptr || length(Gbuffer.emission) > 0)
        return vec3(0.0f); // 暂不计入直接光照

    BRDF brdf(inDir, Gbuffer.shapeNormal, Gbuffer.surfaceNormal);

    float sampleFactor;
    vec3 newDirection = brdf.sample(
            renderData.generator(0),
            renderData.generator(1),
            sampleFactor
    );
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

void render_pixel(const Model &model, const RenderArgs &args, RenderData &renderData, Image &image, int x, int y) {
    const unsigned int
            width = args.width,
            height = args.height,
            spp_direct = args.spp * P_Direct,
            spp_indirect = args.spp - spp_direct;
    const float
            exposure = args.exposure;
    float
            rayX = x - width / 2.0f,
            rayY = y - height / 2.0f;
    vec3 aim = normalize(args.direction + args.accuracy * (rayX * args.right + rayY * args.up));

    RadianceData
            &radiance_i = image.radiance_i[y * width + x],
            &radiance_d = image.radiance_d[y * width + x];
    GbufferData &Gbuffer = image.Gbuffer[y * width + x];

    for (int T = 0; T < spp_direct; T++) {
        vec3 inradiance = sampleDirectLightFromFirstIntersection(Gbuffer, aim, model, renderData);
        inradiance *= exposure;
        radiance_d.radiance += inradiance;
        radiance_d.Var += dot(inradiance, inradiance);
        radiance_d.sampleCount++;
    }

    vec3 intersection = args.position + aim * Gbuffer.depth;
    for (int T = 0; T < spp_indirect; T++) {
        vec3 inradiance = sampleRayFromFirstIntersection(Gbuffer, aim, intersection, model, renderData);
        inradiance *= exposure;
        radiance_i.radiance += inradiance;
        radiance_i.Var += dot(inradiance, inradiance);
        radiance_i.sampleCount++;
    }

    static const float P_skip = 0.5f;
    if (renderData.sobel() < P_skip) {
        renderData.sobel(0);
        renderData.sobel(1);
    }
}

struct PixelPos {
    int x, y;
};

void getRenderSeq(int xL, int xR, int yL, int yR, std::vector<PixelPos> &pixels) {

    if (xL+1 == xR && yL+1 == yR) {
        pixels.push_back({xL, yL});
        return;
    }
    if (xR - xL > yR - yL) {
        int xM = (xL + xR) / 2;
        getRenderSeq(xL, xM, yL, yR, pixels);
        getRenderSeq(xM, xR, yL, yR, pixels);
    }
    else {
        int yM = (yL + yR) / 2;
        getRenderSeq(xL, xR, yL, yM, pixels);
        getRenderSeq(xL, xR, yM, yR, pixels);
    }
}

void render(const Model &model, const RenderArgs &args,
            RenderData &renderData, Image &image, std::atomic<unsigned int> &renderedPixels, PixelPos *pL, PixelPos *pR) {
    static const int PixelsPerReport = 4096;
    static const int PixelsPerAdd = 256;
    for (PixelPos *p = pL; p != pR; p++) {
        render_pixel(model, args, renderData, image, p->x, p->y);
        if ((p-pL) % PixelsPerAdd == 0) {
            renderedPixels += PixelsPerAdd;
            if (renderedPixels % PixelsPerReport == 0)
                std::cout << renderedPixels << "/" << args.width * args.height << "pixels rendered. ("
                          << 1e2f * renderedPixels / (args.width * args.height) << "%)" << std::endl;
        }
    }
}

void render_multiThread(Model &model, const RenderArgs &args) {
    const unsigned int
            startTime = clock(),
            width = args.width,
            height = args.height,
            threads = args.threads;

    std::vector<PixelPos> pixels;
    getRenderSeq(0, width, 0, height, pixels);

    std::cout << "Rendering started with " << threads << " threads." << std::endl;

    Image image(width, height);
    rayCasting(model, args, image);

    std::cout << "Ray casting completed." << std::endl;

    std::atomic<unsigned int> renderedPixels(0);
    std::vector<RenderData> datas(threads);
    std::vector<std::thread> threadPool;

    for (int i = 0; i < threads; i++)
        threadPool.emplace_back([&](int index) {
            render(model, args, datas[index], image, renderedPixels,
                        &pixels.front() + index * pixels.size() / threads,
                        &pixels.front() + (index + 1) * pixels.size() / threads);
        }, i);
    for (auto &thread : threadPool)
        thread.join();

    int T_ray = 0;
    for (auto &data : datas)
        T_ray += data.T_RayAndTexture;

    std::cout << "Rendering completed in " << clock() - startTime << " ms." << std::endl;
    std::cout << "Ray intersection & Texture query time: " << T_ray << " ms*thread." << std::endl;

    image.normalizeRadiance();
    image.shade();
    image.gammaCorrection();
    image.save((args.savePath+"(unfiltered).png").c_str());
    image.filter();
    image.shade();
    image.gammaCorrection();
    image.save((args.savePath+"(filtered).png").c_str());

    std::cout << "Post processing finished. Total: " << clock() - startTime << " ms." << std::endl;
}