#include <thread>
#include <iostream>
#include "render.h"
#include "geometry.h"
#include "image.h"
#include "sampling.h"

RenderData::RenderData(int seed) : gen(seed) {
    C_lightSamples = 0;
}

const float P_RR = 0.5f;
const int maxRayDepth = 8;

vec3 sampleRay(Ray ray, const RayDifferential &base_diff, const Model &model, RenderData &renderData, int &fails, int depth) {
    HitRecord hit = model.rayHit(ray);
    RayDifferential next_diff;
    if (hit.t_max == INFINITY)
        return vec3(0.0f); // No intersection

    BRDF brdf(-ray.direction, ray.origin + ray.direction * hit.t_max); // This is to be initialized by the hit information later

    const auto face = *hit.face;
    const vec3 baryCoords = barycentric(face.v[0], face.v[1], face.v[2], brdf.surface.position);
    const Material& material = *face.material;

    // Initialize the hit information Stage 1, we obtain the normal.
    getHitAllNormals(face, ray.direction, baryCoords, brdf.surface.shapeNormal, brdf.surface.surfaceNormal);

    vec3 hit_dPdx, hit_dPdy; // Derivatives of the position (Calculated outside the if-else block)
    vec3 hit_dDdx, hit_dDdy; // Derivatives of the direction (To be calculated outside the if-else block)

    // Here we just use the normal differential of the shape, which is trivial
    vec3 hit_dNdx = glm::vec3{0.0f}, hit_dNdy = glm::vec3{0.0f}; // Derivatives of the normal

    const auto &hit_t = hit.t_max;
    float dtdx = - dot(base_diff.dPdx + hit_t * base_diff.dDdx, brdf.surface.shapeNormal) / dot(ray.direction, brdf.surface.shapeNormal);
    float dtdy = - dot(base_diff.dPdy + hit_t * base_diff.dDdy, brdf.surface.shapeNormal) / dot(ray.direction, brdf.surface.shapeNormal);
    hit_dPdx = base_diff.dPdx + dtdx * ray.direction + hit_t * hit_dDdx;
    hit_dPdy = base_diff.dPdy + dtdy * ray.direction + hit_t * hit_dDdy;

    // Initialize the hit information Stage 2, we obtain the texture and the tangent space.
    getHitTexture(face, baryCoords, hit_dPdx, hit_dPdy, brdf.surface);
    brdf.genTangentSpace();

    if (glm::length(brdf.surface.emission) > 0)
        return vec3(0.0f); // 不再计算直接光照

    if (depth == maxRayDepth || renderData.gen() > P_RR) {
        renderData.C_lightSamples++;
        vec3 brdfPdf;
        vec3 light = sampleDirectLight(brdf, model, renderData.gen, brdfPdf, fails);
        return light * brdfPdf / ((depth < maxRayDepth) ? (1.0f - P_RR) : 1.0f);
        // RR 停止时，用光源重要性采样计算直接光照
    }

    vec3 newDirection, brdfPdf;
    brdf.sample(renderData.gen, newDirection, brdfPdf, fails); // Sample the new direction
    Ray newRay = {brdf.surface.position, newDirection};

    if (isnan(newDirection))
        return vec3(0.0f); // Invalid direction, then return black

    if (dot(newDirection, brdf.surface.shapeNormal) < 0.0f) {
        // Refraction
        // TODO: After finish the transparent material, we should add the refraction of transparent material

        throw std::runtime_error("Refraction is not supported yet.");

    } else if (dot(newDirection, brdf.surface.shapeNormal) > 0.0f) {
        // Reflection / Diffuse
        // TODO: Calculate the differential of the normal

        // Propagate the differential through the reflection
        float dDNdx = dot(hit_dNdx, ray.direction) + dot(base_diff.dDdx, brdf.surface.shapeNormal);
        float dDNdy = dot(hit_dNdy, ray.direction) + dot(base_diff.dDdy, brdf.surface.shapeNormal);
        hit_dDdx = base_diff.dDdx - 2.0f * (dot(ray.direction, brdf.surface.shapeNormal) * hit_dNdx + dDNdx * brdf.surface.shapeNormal);
        hit_dDdy = base_diff.dDdy - 2.0f * (dot(ray.direction, brdf.surface.shapeNormal) * hit_dNdy + dDNdy * brdf.surface.shapeNormal);

        next_diff = {hit_dPdx, hit_dPdy, hit_dDdx, hit_dDdy};
    }

    vec3 irradiance = sampleRay(newRay, next_diff, model, renderData, fails, depth + 1) / P_RR;

    return brdfPdf * irradiance;
}

void accumulateInwardRadiance(RadianceData &radiance, vec3 inradiance) {
    if (isnan(inradiance)) {
        // inradiance = vec3(0.0f);
        throw std::runtime_error("NaN detected!");
        //std::cerr << "Warning: NaN detected! (accumulate)" << std::endl;
    }
    radiance.radiance += inradiance;
    radiance.Var += dot(inradiance, inradiance);
}

void accumulateInwardRadiance(const BRDF &brdf, const vec3 brdfPdf, const vec3 &light,
                              RadianceData &radiance_d, RadianceData &radiance_s) {

    if (light == vec3(0.0f))
        return;
    vec3 baseColor = (vec3)brdf.surface.baseColor;
    if (length(baseColor) == 0.0f) {
        accumulateInwardRadiance(radiance_s, light * brdfPdf);
        return;
    }
    vec3 baseColor0 = normalize(baseColor);
    static const vec3 White = normalize(vec3(1.0f));
    float XdotY = dot(baseColor0, White);
    if (XdotY > 0.99f) {
        accumulateInwardRadiance(radiance_d, light * brdfPdf / baseColor);
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
    accumulateInwardRadiance(radiance_d, light * B / length(baseColor));
    accumulateInwardRadiance(radiance_s, light * (brdfPdf - brdfPdf_base));
}

void sampleDirectLightFromFirstIntersection(const HitInfo &hitInfo, const vec3 &inDir, const Model &model,
                                            RenderData &renderData, RadianceData &radiance_Dd, RadianceData &radiance_Ds, int &fails) {

    if (isnan(hitInfo.position) || length(hitInfo.emission) > 0)
        return; // 暂不计入直接光照

    BRDF brdf(-inDir, hitInfo);
    brdf.genTangentSpace();
    vec3 brdfPdf;
    vec3 light = sampleDirectLight(brdf, model, renderData.gen, brdfPdf, fails);
    accumulateInwardRadiance(brdf, brdfPdf, light, radiance_Dd, radiance_Ds);
}

void sampleRayFromFirstIntersection(const HitInfo &hitInfo, const vec3 &inDir, const RayDifferential base_diff, const Model &model,
                                    RenderData &renderData, RadianceData &radiance_Id, RadianceData &radiance_Is, int &fails) {

    if (isnan(hitInfo.position) || length(hitInfo.emission) > 0)
        return; // 暂不计入直接光照

    BRDF brdf(-inDir, hitInfo);
    brdf.genTangentSpace();
    vec3 newDirection, brdfPdf;
    brdf.sample(renderData.gen, newDirection, brdfPdf, fails);
    if (isnan(newDirection))
        return;
    Ray ray = {hitInfo.position, newDirection};
    vec3 light = sampleRay(ray, base_diff, model, renderData, fails, 1);
    accumulateInwardRadiance(brdf, brdfPdf, light, radiance_Id, radiance_Is);
}

void rayCasting(Model &model, const RenderArgs &args, Image &image) {
    const int
            width = args.width,
            height = args.height;
    vec3
            view = args.direction,
            right = args.right,
            up = args.up,
            position = args.position;
    const float
            accuracy = args.accuracy;
    for (int x = 0; x < width; x++)
        for (int y = 0; y < height; y++) {

            float
                    rayX = static_cast<float>(x) - width / 2.0f,
                    rayY = static_cast<float>(y) - height / 2.0f;
            vec3 d = view + accuracy * (rayX * right + rayY * up);

            vec3 D = normalize(d);
            Ray ray = {position, D};

            vec3 base_dPdx = vec3(0.0f), base_dPdy = vec3(0.0f);
            vec3 base_dDdx, base_dDdy;

            vec3 dddx = accuracy * right;
            vec3 dddy = accuracy * up;

            float d_dot_d = dot(d, d);
            float d_dot_dddx = dot(d, dddx);
            float d_dot_dddy = dot(d, dddy);

            base_dDdx = (d_dot_d * dddx - d_dot_dddx * d) / (sqrt(d_dot_d) * d_dot_d);
            base_dDdy = (d_dot_d * dddy - d_dot_dddy * d) / (sqrt(d_dot_d) * d_dot_d);

            HitRecord hit = model.rayHit(ray);
            HitInfo &Gbuffer = image.Gbuffer[y * width + x];

            if (hit.face == nullptr)
                continue;

            Gbuffer.position = ray.origin + ray.direction * hit.t_max;

            const vec3 baryCoords = barycentric(hit.face->v[0], hit.face->v[1], hit.face->v[2], Gbuffer.position);
            getHitAllNormals(*hit.face, ray.direction, baryCoords, Gbuffer.shapeNormal, Gbuffer.surfaceNormal);

            vec3 hit_dPdx, hit_dPdy; // Derivatives of the position, a trivial solution

            const auto &hit_t = hit.t_max;
            float dtdx = - dot(hit_t * base_dDdx + base_dPdx, Gbuffer.shapeNormal) / dot(ray.direction, Gbuffer.shapeNormal),
                  dtdy = - dot(hit_t * base_dDdy + base_dPdy, Gbuffer.shapeNormal) / dot(ray.direction, Gbuffer.shapeNormal);
            hit_dPdx = base_dPdx + dtdx * ray.direction + hit_t * base_dDdx;
            hit_dPdy = base_dPdy + dtdy * ray.direction + hit_t * base_dDdy;

            getHitTexture(*hit.face, baryCoords, hit_dPdx, hit_dPdy, Gbuffer);
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
    HitInfo &Gbuffer = image.Gbuffer[id];
    if (isnan(Gbuffer.position))
        return;
    const int
            spp_direct = args.spp * args.P_Direct,
            spp_indirect = args.spp - spp_direct;
    const float
            exposure = args.exposure;
    vec3 d = Gbuffer.position - args.position;
    vec3 D = normalize(d); // D
    RadianceData
            &radiance_Dd = image.radiance_Dd[id],
            &radiance_Ds = image.radiance_Ds[id],
            &radiance_Id = image.radiance_Id[id],
            &radiance_Is = image.radiance_Is[id];
    int trys_direct = 0;

    RayDifferential base_diff;

    base_diff.dPdx = vec3(0.0f), base_diff.dPdy = vec3(0.0f);
    const vec3 dddx = args.accuracy * args.right, dddy = args.accuracy * args.up;
    float d_dot_d = dot(d, d), d_dot_dddx = dot(d, dddx), d_dot_dddy = dot(d, dddy);

    base_diff.dDdx = (d_dot_d * dddx - d * d_dot_dddx) / (sqrt(d_dot_d) * d_dot_d);
    base_diff.dDdy = (d_dot_d * dddy - d * d_dot_dddy) / (sqrt(d_dot_d) * d_dot_d);

    for (int T = 0; T < spp_direct; T++)
        sampleDirectLightFromFirstIntersection(Gbuffer, D, model,
                                               renderData, radiance_Dd, radiance_Ds, trys_direct);
    trys_direct += spp_direct;
    normalize(radiance_Dd, exposure, trys_direct);
    normalize(radiance_Ds, exposure, trys_direct);
    int trys_indirect = 0;
    for (int T = 0; T < spp_indirect; T++)
        sampleRayFromFirstIntersection(Gbuffer, D, base_diff, model,
                                       renderData, radiance_Id, radiance_Is, trys_indirect);
    trys_indirect += spp_indirect;
    normalize(radiance_Id, exposure, trys_indirect);
    normalize(radiance_Is, exposure, trys_indirect);
}

/*
 * Render Workers
 * Master call this worker function to render the scene.
 */
void renderWorker(const Model &model, const RenderArgs &args,
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

    auto exportImage = [&](const std::string &tag, int shadeOptions) {
        image.shade(exposure, shadeOptions);
        image.save((args.savePath + "(" + tag + ").png").c_str());
    };

    exportImage("DiffuseColor", Image::BaseColor);
    exportImage("shapeNormal", Image::shapeNormal);
    exportImage("surfaceNormal", Image::surfaceNormal);

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
            renderWorker(model, args, datas[index], image,
                         index * width / threads, (index + 1) * width / threads);
        }, i);
    for (auto &thread : threadPool)
        thread.join();

    int C_lightSamples = 0;
    for (const auto &data : datas)
        C_lightSamples += data.C_lightSamples;

    std::cout << "Rendering completed in " << clock() - startTime << " ms." << std::endl;
    std::cout << "Direct light samples: " << C_lightSamples << std::endl;

    exportImage("Direct_Diffuse", Image::Direct_Diffuse);
    exportImage("Direct_Specular", Image::Direct_Specular);
    exportImage("Indirect_Diffuse", Image::Indirect_Diffuse);
    exportImage("Indirect_Specular", Image::Indirect_Specular);
    exportImage("Raw", Image::Full);
    image.filter();
    exportImage("Direct_Diffuse_Filter", Image::Direct_Diffuse);
    exportImage("Direct_Specular_Filter", Image::Direct_Specular);
    exportImage("Indirect_Diffuse_Filter", Image::Indirect_Diffuse);
    exportImage("Indirect_Specular_Filter", Image::Indirect_Specular);
    exportImage("Filter", Image::Full);

    std::cout << "Post processing finished. Total: " << clock() - startTime << " ms." << std::endl;
}