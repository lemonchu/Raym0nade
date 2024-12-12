#include <thread>
#include <random>
#include <iostream>
#include "render.h"
#include "image.h"
#include "sampling.h"

RenderData::RenderData(int seed) : gen(seed) {
    T_RayAndTexture = C_lightSamples = 0;
}

const float P_RR = 0.5f;
const int maxRayDepth = 8;

vec3 sampleRay(Ray ray, const Model &model, RenderData &renderData, int &fails, int depth) {

    renderData.T_RayAndTexture -= clock();
    BRDF brdf(-ray.direction);
    model.rayHit(ray, brdf.surface);
    renderData.T_RayAndTexture += clock();

    if (brdf.surface.t == INFINITY)
        return vec3(0.0f);

    if (length(brdf.surface.emission) > 0)
        return vec3(0.0f); // 不再计算直接光照

    vec3 intersection = ray.origin + ray.direction * brdf.surface.t;
    if (depth == maxRayDepth || renderData.gen() > P_RR) {
        renderData.C_lightSamples++;
        vec3 brdfPdf;
        vec3 light = sampleDirectLight(intersection, brdf, model, renderData.gen, brdfPdf, fails);
        return light * brdfPdf / ((depth < maxRayDepth) ? (1.0f - P_RR) : 1.0f);
        // RR 停止时，用光源重要性采样计算直接光照
    }

    brdf.genTangentSpace();
    vec3 newDirection, brdfPdf;
    brdf.sample(renderData.gen, newDirection, brdfPdf, fails);
    if (isnan(newDirection))
        return vec3(0.0f);
    ray = {intersection, newDirection};
    vec3 inradiance = sampleRay(ray, model, renderData, fails, depth + 1) / P_RR;
    return brdfPdf * inradiance;
}

void accumulateInradiance(RadianceData &radiance, vec3 inradiance) {
    if (isnan(inradiance)) {
        inradiance = vec3(0.0f);
        std::cerr << "Warning: NaN detected! (accumulate)" << std::endl;
    }
    radiance.radiance += inradiance;
    radiance.Var += dot(inradiance, inradiance);
}

void accumulateInradiance(const BRDF &brdf, const vec3 brdfPdf, const vec3 &light,
                          RadianceData &radiance_d, RadianceData &radiance_s) {

    if (light == vec3(0.0f))
        return;
    vec3 baseColor = (vec3)brdf.surface.baseColor;
    if (length(baseColor) == 0.0f) {
        accumulateInradiance(radiance_s, light * brdfPdf);
        return;
    }
    vec3 baseColor0 = normalize(baseColor);
    static const vec3 White = normalize(vec3(1.0f));
    float XdotY = dot(baseColor0, White);
    if (XdotY > 0.99f) {
        accumulateInradiance(radiance_d, light * brdfPdf / baseColor);
        return ;
    }
    // Consider Pdf = a * white + b * baseColor
    vec3 perp = normalize(cross(baseColor0, White));
    vec3 brdfPdf_p = brdfPdf - perp * dot(perp, brdfPdf);
    float d1 = dot(brdfPdf_p, White);
    float d2 = dot(brdfPdf_p, baseColor0);
    float AplusB = (d1+d2) / (1+XdotY);
    float AminusB = (d1-d2) / (1-XdotY);
    float B = (AplusB - AminusB) / 2.0f;
    vec3 brdfPdf_base = B * baseColor0;
    accumulateInradiance(radiance_d, light * B / length(baseColor));
    accumulateInradiance(radiance_s, light * (brdfPdf - brdfPdf_base));
}

void sampleDirectLightFromFirstIntersection(const GbufferData &Gbuffer, const vec3 &inDir, const vec3 &emission, const Model &model,
                                            RenderData &renderData, RadianceData &radiance_Dd, RadianceData &radiance_Ds, int &fails) {

    if (Gbuffer.face == nullptr || length(emission) > 0)
        return; // 暂不计入直接光照

    BRDF brdf(-inDir);
    getHitInfo(*Gbuffer.face, Gbuffer.position, inDir, brdf.surface);
    brdf.genTangentSpace();
    vec3 brdfPdf;
    vec3 light = sampleDirectLight(Gbuffer.position, brdf, model, renderData.gen, brdfPdf, fails);
    accumulateInradiance(brdf, brdfPdf, light, radiance_Dd, radiance_Ds);
}


void sampleRayFromFirstIntersection(const GbufferData &Gbuffer, const vec3 &inDir, const vec3 &emission, const Model &model,
                                    RenderData &renderData, RadianceData &radiance_Id, RadianceData &radiance_Is, int &fails) {

    if (Gbuffer.face == nullptr || length(emission) > 0)
        return; // 暂不计入直接光照

    BRDF brdf(-inDir);
    getHitInfo(*Gbuffer.face, Gbuffer.position, inDir, brdf.surface);
    brdf.genTangentSpace();
    vec3 newDirection, brdfPdf;
    brdf.sample(renderData.gen, newDirection, brdfPdf, fails);
    if (isnan(newDirection))
        return;
    Ray ray = {Gbuffer.position, newDirection};
    vec3 light = sampleRay(ray, model, renderData, fails, 1);
    accumulateInradiance(brdf, brdfPdf, light, radiance_Id, radiance_Is);
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

void normalize(RadianceData &radiance, float exposure, int spp) {
    radiance.radiance *= exposure;
    radiance.Var *= exposure * exposure;
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
    RadianceData
            &radiance_Dd = image.radiance_Dd[id],
            &radiance_Ds = image.radiance_Ds[id],
            &radiance_Id = image.radiance_Id[id],
            &radiance_Is = image.radiance_Is[id];
    HitInfo hitInfo;
    getHitInfo(*Gbuffer.face, Gbuffer.position, aim, hitInfo);
    int trys_direct = 0;
    for (int T = 0; T < spp_direct; T++)
        sampleDirectLightFromFirstIntersection(Gbuffer, aim, hitInfo.emission, model,
                                               renderData, radiance_Dd, radiance_Ds, trys_direct);
    trys_direct += spp_direct;
    normalize(radiance_Dd, exposure, trys_direct);
    normalize(radiance_Ds, exposure, trys_direct);
    int trys_indirect = 0;
    for (int T = 0; T < spp_indirect; T++)
        sampleRayFromFirstIntersection(Gbuffer, aim, hitInfo.emission, model,
                                       renderData, radiance_Id, radiance_Is, trys_indirect);
    trys_indirect += spp_indirect;
    normalize(radiance_Id, exposure, trys_indirect);
    normalize(radiance_Is, exposure, trys_indirect);
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
    const vec3 position = args.position;
    const float exposure = args.exposure;

    std::cout << "Rendering started with " << threads << " threads." << std::endl;

    Image image(width, height);
    rayCasting(model, args, image);

    std::cout << "Ray casting completed. (" << clock()-startTime << " ms)"<< std::endl;

    image.shade(position, exposure, Image::BaseColor);
    image.save((args.savePath+"(DiffuseColor).png").c_str());
    image.shade(position, exposure, Image::shapeNormal);
    image.save((args.savePath+"(shapeNormal).png").c_str());
    image.shade(position, exposure, Image::surfaceNormal);
    image.save((args.savePath+"(surfaceNormal).png").c_str());

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


    image.shade(position, exposure, Image::DirectLight | Image::Diffuse);
    image.save((args.savePath+"(Direct_Diffuse).png").c_str());
    image.shade(position, exposure, Image::DirectLight | Image::Specular);
    image.save((args.savePath+"(Direct_Specular).png").c_str());
    image.shade(position, exposure, Image::IndirectLight | Image::Diffuse);
    image.save((args.savePath+"(Indirect_Diffuse).png").c_str());
    image.shade(position, exposure, Image::IndirectLight | Image::Specular);
    image.save((args.savePath+"(Indirect_Specular).png").c_str());
    image.shade(position, exposure,
                Image::DirectLight | Image::IndirectLight |
                Image::Diffuse | Image::Specular |
                Image::Emission | Image::BaseColor);
    image.save((args.savePath+"(raw).png").c_str());
    image.filter();
    image.shade(position, exposure, Image::DirectLight | Image::Diffuse);
    image.save((args.savePath+"(Direct_Diffuse_Filter).png").c_str());
    image.shade(position, exposure, Image::DirectLight | Image::Specular);
    image.save((args.savePath+"(Direct_Specular_Filter).png").c_str());
    image.shade(position, exposure, Image::IndirectLight | Image::Diffuse);
    image.save((args.savePath+"(Indirect_Diffuse_Filter).png").c_str());
    image.shade(position, exposure, Image::IndirectLight | Image::Specular);
    image.save((args.savePath+"(Indirect_Specular_Filter).png").c_str());
    image.shade(position, exposure,
                Image::DirectLight | Image::IndirectLight |
                Image::Diffuse | Image::Specular |
                Image::Emission | Image::BaseColor);
    image.save((args.savePath+"(Filter).png").c_str());

    std::cout << "Post processing finished. Total: " << clock() - startTime << " ms." << std::endl;
}