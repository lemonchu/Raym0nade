#include <mutex>
#include <thread>
#include <iostream>
#include "render.h"
#include "geometry.h"
#include "image.h"
#include "sampling.h"

RenderData::RenderData(int seed) : gen(seed) {
    C_lightSamples = 0;
}

MediumData::MediumData(float ior, const vec3 &absorb) : ior(ior), absorb(absorb) {}

void Medium::Init() {
    mediums.clear();
    mediums.emplace(-1, MediumData(1.0f, vec3(1.0f)));
}

void Medium::insert(int id, float ior, const vec3 &absorb) {
    mediums.emplace(id, MediumData(ior, absorb));
}

void Medium::erase(int id) {
    mediums.erase(id);
}

float Medium::ior() const {
    float ior = 1.0f;
    for (const auto &medium: mediums)
        ior = std::max(ior, medium.second.ior);
    return ior;
}

vec3 Medium::absorb() const {
    vec3 absorb = vec3(1.0f);
    for (const auto &medium: mediums)
        absorb *= medium.second.absorb;
    return absorb;
}

int Medium::size() const {
    return static_cast<int>(mediums.size());
}

void calc_dPdxy(const Ray &ray, float hit_t, const vec3 &normal, const RayDifferential &base_diff,
                vec3 &hit_dPdx, vec3 &hit_dPdy) {
    float dtdx = -dot(base_diff.dPdx + hit_t * base_diff.dDdx, normal) / dot(ray.direction, normal);
    float dtdy = -dot(base_diff.dPdy + hit_t * base_diff.dDdy, normal) / dot(ray.direction, normal);
    hit_dPdx = base_diff.dPdx + dtdx * ray.direction + hit_t * base_diff.dDdx;
    hit_dPdy = base_diff.dPdy + dtdy * ray.direction + hit_t * base_diff.dDdy;
}

const vec3 hit_dNdx = glm::vec3{0.0f}, hit_dNdy = glm::vec3{0.0f}; // Derivatives of the normal
void calc_dDdxy(const Ray &ray, const vec3 &normal, const RayDifferential &base_diff,
                vec3 &hit_dDdx, vec3 &hit_dDdy) {

    float dDNdx = dot(hit_dNdx, ray.direction) + dot(base_diff.dDdx, normal);
    float dDNdy = dot(hit_dNdy, ray.direction) + dot(base_diff.dDdy, normal);
    hit_dDdx = base_diff.dDdx - 2.0f * (dot(ray.direction, normal) * hit_dNdx + dDNdx * normal);
    hit_dDdy = base_diff.dDdy - 2.0f * (dot(ray.direction, normal) * hit_dNdy + dDNdy * normal);
}

void getHitInfo(const HitRecord &hit, const Ray &ray, const RayDifferential &base_diff,
                vec3 &hit_dPdx, vec3 &hit_dPdy, HitInfo &hitInfo) {
    const auto face = *hit.face;
    const vec3 baryCoords = barycentric(face.v[0], face.v[1], face.v[2], hitInfo.position);

    // Initialize the hit information Stage 1, we obtain the raw normal.
    getHitNormals(face, ray.direction, baryCoords,
                  hitInfo.shapeNormal, hitInfo.surfaceNormal, hitInfo.entering);
    vec3 surfaceNormal_raw = hitInfo.surfaceNormal;

    // Initialize the hit information Stage 2, we obtain the texture and the tangent space.
    calc_dPdxy(ray, hit.t_max, hitInfo.shapeNormal, base_diff, hit_dPdx, hit_dPdy);
    getHitMaterial(face, baryCoords, hit_dPdx, hit_dPdy, hitInfo);
    if (dot(hitInfo.surfaceNormal, ray.direction) >= 0.0f)
        hitInfo.surfaceNormal = surfaceNormal_raw;
    if (dot(hitInfo.surfaceNormal, ray.direction) >= 0.0f)
        hitInfo.surfaceNormal = hitInfo.shapeNormal;
}

// #define DEBUG_sampleRay

vec3 getAbsorb(const vec3 &absorb, float dis) { // 计算透射颜色
    static const float omega_absorb = 32.0f;
    float Clum = dot(absorb, RGB_Weight);
    return (Clum < 1.0f - eps_zero) ? pow_s(absorb, omega_absorb * dis) : vec3(1.0f);
}

void calcEta(Medium &mediums, HitInfo &hitInfo) { // 计算相对折射率
    float eta1 = mediums.ior(), eta2;
    if (hitInfo.entering)
        eta2 = std::max(eta1, hitInfo.eta);
    else {
        mediums.erase(hitInfo.id);
        eta2 = mediums.ior();
        mediums.insert(hitInfo.id, hitInfo.eta, hitInfo.baseColor);
    }
    hitInfo.eta = eta1 / eta2;
}

void scaling(std::vector<LightSample> &samples, float factor) {
    for (auto &sample: samples)
        sample.bsdfPdf *= factor;
}

void passBsdf(std::vector<LightSample> &samples, const vec3 &bsdfPdf) {
    for (auto &sample: samples) {
        sample.light *= sample.bsdfPdf;
        sample.bsdfPdf = bsdfPdf;
    }
}

void mulWeight(std::vector<LightSample> &samples, float weight) {
    for (auto &sample: samples)
        sample.weight *= weight;
}

static float P0_reflect = 0.24f;
const float regularizationFactor = 1.0f;

std::vector<LightSample> sampleRay
        (Ray ray, const RayDifferential &base_diff, const Model &model,
         RenderData &renderData, float roughnessFactor, bool excludeDirectLight, int depth) {

    static const int maxRayDepth = 16;
    static const int sampleCount[maxRayDepth + 1] = {0, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6};

    HitRecord hit = model.rayHit(ray);
    if (hit.t_max == INFINITY) {
        if (excludeDirectLight || model.skyMap.empty())
            return {};
        return {LightSample{vec3(1.0f), model.skyMap.get(ray.direction), 1.0f}};
    }

    // 获取表面信息
    BSDF bsdf(-ray.direction, ray.origin + ray.direction * hit.t_max);
    vec3 hit_dPdx, hit_dPdy;
    getHitInfo(hit, ray, base_diff, hit_dPdx, hit_dPdy, bsdf.surface);

    float ior = bsdf.surface.eta; // 保存绝对折射率

    // 路径正则化
    roughnessFactor = std::max(roughnessFactor, regularizationFactor * bsdf.surface.roughness);
    bsdf.surface.roughness = std::max(bsdf.surface.roughness, roughnessFactor);

    // 根据光滑程度决定 RR 的截止概率
    static constexpr float P_RR_Rough = 0.5f, P_RR_Smooth = 1.0f;
    const float P_RR = P_RR_Smooth + (P_RR_Rough - P_RR_Smooth) * sqrt(bsdf.surface.roughness);
    static constexpr float P_RR_Therehold = 0.9f;
    bool doDirectLightSample = P_RR < P_RR_Therehold;

    // 计算介质的吸收率
    vec3 absorb = getAbsorb(renderData.mediums.absorb(), hit.t_max);

#ifdef DEBUG_sampleRay
    std::cout << "depth: " << depth << std::endl;
    std::cout << "entering: " << (bsdf.surface.entering ? "true" : "false") << std::endl ;
    std::cout << "name: " << hit.face->material->name << std::endl;
#endif

    // 不计直接光照
    if (length(bsdf.surface.emission) > eps_zero) {
        if (excludeDirectLight)
            return {};
        return {LightSample{absorb, bsdf.surface.emission, 1.0f}};
    }

    float P_reflect = 1.0f, F;
    vec3 refractDir(NAN);

    // 透明材质，计算反射和折射的比率（顺便计算折射方向）
    if (bsdf.surface.opacity < eps_zero) {
        calcEta(renderData.mediums, bsdf.surface);
        bsdf.preciseRefraction(refractDir, F);
        if (renderData.mediums.size() == 1) // 在空气中，增加反射采样率
            P_reflect = P0_reflect + (1.0f - P0_reflect) * F;
        else
            P_reflect = F;
        P_reflect = std::max(P_reflect - 1e-3f, 0.0f);
#ifdef DEBUG_sampleRay
        std::cout << "    eta: " << bsdf.surface.eta << std::endl;
        std::cout << "    refractDir: " << refractDir.x << " " << refractDir.y << " " << refractDir.z << std::endl;
        std::cout << "    P_reflect: " << P_reflect << std::endl;
#endif
    }

#ifdef DEBUG_sampleRay
    std::cout << "        Absorb_depth: " << hit.t_max << std::endl;
    std::cout << "        Absorb: " << absorb.x << " " << absorb.y << " " << absorb.z << std::endl;
#endif

    vec3 newDir, bsdfPdf(NAN);
    RayDifferential next_diff;
    int fails = 0;

    if (renderData.gen() <= P_reflect) {
        // 发生反射
#ifdef DEBUG_sampleRay
        std::cout << "reflect" << std::endl;
#endif
        doDirectLightSample &=
                bsdf.surface.entering &&
                renderData.mediums.size() == 1;
        // 仅在空气中进行直接光照采样（在空气中反射）

        if (doDirectLightSample) {
            if (depth >= maxRayDepth || renderData.gen() > P_RR) {
                renderData.C_lightSamples++;
                std::vector<LightSample> samples =
                        sampleDirectLight(bsdf, model, renderData.gen, sampleCount[depth]);
                float factor = 1.0f / P_reflect;
                if (depth < maxRayDepth)
                    factor /= (1.0f - P_RR);
                scaling(samples, factor);
#ifdef DEBUG_sampleRay
                std::cout << "bsdfPdf : " << bsdfPdf.x << " " << bsdfPdf.y << " " << bsdfPdf.z << std::endl;
                std::cout << "depth :" << depth << std::endl;
                std::cout << std::endl;
#endif
                return samples;
            }
        }

        // 采样反射方向并获得 brdfPdf
        bsdf.sampleReflection(renderData.gen, newDir, bsdfPdf, fails);
        if (!isfinite(bsdfPdf)) {
            std::cout << "bsdfPdf is NaN! (reflection)" << std::endl;
        }
        bsdfPdf /= P_reflect;
        if (doDirectLightSample)
            bsdfPdf /= P_RR;

        // 计算光微分
        vec3 hit_dDdx, hit_dDdy;
        calc_dDdxy(ray, bsdf.surface.surfaceNormal, base_diff, hit_dDdx, hit_dDdy);
        next_diff = {hit_dPdx, hit_dPdy, hit_dDdx, hit_dDdy};

    } else {
        // 发生折射
#ifdef DEBUG_sampleRay
        std::cout << "refract" << std::endl;
#endif
        // Todo: 折射的 RayDifferential 计算
        doDirectLightSample &=
                !bsdf.surface.entering &&
                renderData.mediums.size() == 2;
        // 仅在空气中进行直接光照采样（折射离开介质，进入空气）

        if (doDirectLightSample) {
            if (depth >= maxRayDepth || renderData.gen() > P_RR) {
                renderData.C_lightSamples++;
                std::vector<LightSample> samples =
                        sampleDirectLight(bsdf, model, renderData.gen, sampleCount[depth]);
                float factor = 1.0f / (1.0f - P_reflect);
                if (depth < maxRayDepth)
                    factor /= (1.0f - P_RR);
                scaling(samples, factor);
                passBsdf(samples, absorb);
                return samples;
            }
        }

        // 取折射方向
        if (abs(bsdf.surface.eta - 1.0f) < eps_zero) {
            // 介质内部，不反射
            newDir = ray.direction;
            bsdfPdf = vec3(1.0f);
            next_diff = base_diff;
        } else {
            bsdf.sampleBTDF(renderData.gen, newDir, bsdfPdf, fails);
            bsdfPdf *= (1.0f - F);
            next_diff = base_diff;
            // Todo: 折射的 RayDifferential 计算
            next_diff = base_diff;
        }

        bsdfPdf /= (1.0f - P_reflect);
        if (doDirectLightSample)
            bsdfPdf /= P_RR;

        if (bsdf.surface.entering)
            renderData.mediums.insert(bsdf.surface.id, ior, bsdf.surface.baseColor);
        else
            renderData.mediums.erase(bsdf.surface.id);
    }

    // 若方向无效或递归过深，返回 0
    if (!isfinite(newDir) || depth == maxRayDepth)
        return {}; // Invalid direction, then return black

    Ray newRay = {bsdf.surface.position, newDir};

    bsdfPdf *= absorb;

#ifdef DEBUG_sampleRay
    std::cout << "newPos: " << newRay.origin.x << " " << newRay.origin.y << " " << newRay.origin.z << std::endl;
    std::cout << "newDir: " << newDir.x << " " << newDir.y << " " << newDir.z << std::endl;
    std::cout << "bsdfPdf: " << bsdfPdf.x << " " << bsdfPdf.y << " " << bsdfPdf.z << std::endl << std::endl;
#endif

    std::vector<LightSample> samples =
            sampleRay(newRay, next_diff, model,
                      renderData, roughnessFactor, doDirectLightSample, depth + 1);

    if (!isfinite(bsdfPdf)) {
        std::cout << "bsdfPdf is NaN! " << depth << std::endl;
    }
    passBsdf(samples, bsdfPdf);
    if (fails > 0)
        mulWeight(samples, 1.0f / static_cast<float>(1 + fails));

    return samples;
}

void sampleIndirectLightFromFirstIntersection(
        const HitInfo &hitInfo, const vec3 &origin, const RayDifferential &base_diff, const Model &model,
        RenderData &renderData, std::vector<LightSample> &samples) {

    // 不计入直接光照
    if (length(hitInfo.emission) > 0)
        return;

    // 初始化介质信息
    renderData.mediums.Init();

    vec3 inDir = normalize(hitInfo.position - origin);
    Ray ray = {origin, inDir};
    float hit_t = length(hitInfo.position - origin);
    BSDF bsdf(-inDir, hitInfo);

    float roughnessFactor = hitInfo.roughness * regularizationFactor;

    float ior = bsdf.surface.eta;
    bsdf.surface.eta = 1.0f / bsdf.surface.eta; // 初始时必然从空气进入介质

    float P_reflect = 1.0f, F;
    vec3 refractDir;
    if (bsdf.surface.opacity < eps_zero) {
        // 透明材质
        bsdf.preciseRefraction(refractDir, F);
        P_reflect = P0_reflect + (1.0f - P0_reflect) * F;
#ifdef DEBUG_sampleRay
        std::cout << "    refractDir: " << refractDir.x << " " << refractDir.y << " " << refractDir.z << std::endl;
        std::cout << "    P_reflect: " << P_reflect << std::endl;
#endif
    }

    vec3 newDir, bsdfPdf = vec3(NAN);
    RayDifferential next_diff;
    int fails = 0;

    if (renderData.gen() < P_reflect) {
        // 反射
#ifdef DEBUG_sampleRay
        std::cout << "reflect" << std::endl;
#endif
        vec3 hit_dPdx, hit_dPdy;
        calc_dPdxy(ray, hit_t, bsdf.surface.shapeNormal, base_diff, hit_dPdx, hit_dPdy);
        vec3 hit_dDdx, hit_dDdy;
        calc_dDdxy(ray, bsdf.surface.surfaceNormal, base_diff, hit_dDdx, hit_dDdy);
        next_diff = {hit_dPdx, hit_dPdy, hit_dDdx, hit_dDdy};

        bsdf.sampleReflection(renderData.gen, newDir, bsdfPdf, fails);
        bsdfPdf /= P_reflect;
        if (isinf(bsdfPdf.x)) {
            std::cerr << "bsdfPdf is INF! (reflection)" << std::endl;
        }
    } else {
        // 折射
        // Todo: 折射的 RayDifferential 计算
#ifdef DEBUG_sampleRay
        std::cout << "refract" << std::endl;
#endif
        bsdf.sampleBTDF(renderData.gen, newDir, bsdfPdf, fails);
        bsdfPdf *= (1.0f - F);
        bsdfPdf /= (1.0f - P_reflect);

        next_diff = base_diff;

        if (bsdf.surface.entering)
            renderData.mediums.insert(bsdf.surface.id, ior, bsdf.surface.baseColor);
        else
            renderData.mediums.erase(bsdf.surface.id);
    }

    if (!isfinite(newDir))
        return;

    Ray newRay = {bsdf.surface.position, newDir};

#ifdef DEBUG_sampleRay
    std::cout << "depth: 0" << std::endl;
    std::cout << "newPos: " << newRay.origin.x << " " << newRay.origin.y << " " << newRay.origin.z << std::endl;
    std::cout << "newDir: " << newDir.x << " " << newDir.y << " " << newDir.z << std::endl;
    std::cout << "bsdfPdf: " << bsdfPdf.x << " " << bsdfPdf.y << " " << bsdfPdf.z << std::endl;
    std::cout << "entering: " << (bsdf.surface.entering ? "true" : "false") << std::endl << std::endl;
#endif

    std::vector<LightSample> samplesNew =
            sampleRay(newRay, next_diff, model,
                      renderData, roughnessFactor, true, 1);

    passBsdf(samplesNew, bsdfPdf);
    if (fails > 0)
        mulWeight(samplesNew, 1.0f / static_cast<float>(1 + fails));

#ifdef DEBUG_sampleRay
    for (const auto &sample : samplesNew) {
        if (length(sample.light) > 0.1)
            std::cout << "!!!" << std::endl;
        std::cout << "    light: " << sample.light.x << " " << sample.light.y << " " << sample.light.z << std::endl;
        std::cout << "    bsdfPdf: " << sample.bsdfPdf.x << " " << sample.bsdfPdf.y << " " << sample.bsdfPdf.z << std::endl;
        std::cout << "    weight: " << sample.weight << std::endl;
    }
    std::cout << std::endl << std::endl;
#endif

    if (!isfinite(bsdfPdf)) {
        std::cerr << "bsdfPdf is NAN! (final)" << std::endl;
    }

    for (const auto &sample: samplesNew)
        samples.push_back(sample);
}

void sampleDirectLightFromFirstIntersection(const HitInfo &hitInfo, const vec3 &origin, const Model &model, int spp,
                                            RenderData &renderData, RadianceData &radiance_Dd,
                                            RadianceData &radiance_Ds) {

    if (!isfinite(hitInfo.position) || length(hitInfo.emission) > 0)
        return; // 暂不计入直接光照

    vec3 inDir = normalize(hitInfo.position - origin);
    BSDF bsdf(-inDir, hitInfo);
    std::vector<LightSample> samples =
            sampleDirectLight(bsdf, model, renderData.gen, spp);
    for (const auto &sample: samples)
        accumulateInwardRadiance(hitInfo.baseColor, sample, radiance_Dd, radiance_Ds);
}

void initRayDiff(const vec3 &d, const RenderArgs &args, RayDifferential &base_diff) {
    base_diff.dPdx = vec3(0.0f), base_diff.dPdy = vec3(0.0f);
    const vec3 dddx = args.accuracy * args.right, dddy = args.accuracy * args.up;
    float d_dot_d = dot(d, d), d_dot_dddx = dot(d, dddx), d_dot_dddy = dot(d, dddy);

    base_diff.dDdx = (d_dot_d * dddx - d * d_dot_dddx) / (sqrt(d_dot_d) * d_dot_d);
    base_diff.dDdy = (d_dot_d * dddy - d * d_dot_dddy) / (sqrt(d_dot_d) * d_dot_d);
}

void renderPixel(const Model &model, const RenderArgs &args,
                 RenderData &renderData, Image &image, int x, int y) {
    //renderData.gen.mt.seed(x^y);
#ifdef DEBUG_sampleRay
    if (x!=373||y!=30)
        return ;
#endif
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

    int id = y * width + x;
    float
            rayX = float(x) - float(width) / 2.0f,
            rayY = float(y) - float(height) / 2.0f;
    vec3 d = view + accuracy * (rayX * right + rayY * up);
    HitInfo &Gbuffer = image.Gbuffer[id];

    RayDifferential base_diff;
    initRayDiff(d, args, base_diff);

    vec3 Dir = normalize(d);
    Ray ray = {position, Dir};

    HitRecord hit = model.rayHit(ray);
    if (hit.t_max == INFINITY) {
        Gbuffer.position = vec3(NAN);
        if (!model.skyMap.empty())
            Gbuffer.emission = model.skyMap.get(Dir);
        return;
    }

    Gbuffer.position = ray.origin + hit.t_max * ray.direction;
    vec3 hit_dPdx, hit_dPdy;
    getHitInfo(hit, ray, base_diff, hit_dPdx, hit_dPdy, Gbuffer);
    bool opacity = (Gbuffer.opacity > 1.0f - eps_zero);
    static const int MultiSampleForRefraction = 16;
    const int
            spp_direct = args.spp * args.P_Direct,
            spp_indirect = (args.spp - spp_direct) * (opacity ? 1 : MultiSampleForRefraction);
    const float
            exposure = args.exposure;
    RadianceData
            &radiance_Dd = image.radiance_Dd[id],
            &radiance_Ds = image.radiance_Ds[id],
            &radiance_Id = image.radiance_Id[id],
            &radiance_Is = image.radiance_Is[id];

    auto calcVar = [](RadianceData &radiance, float exposure) -> void {
        radiance.radiance *= exposure;
        radiance.Var *= exposure * exposure;
        radiance.Var -= dot(radiance.radiance, radiance.radiance);
        if (radiance.Var < 0.0f)
            radiance.Var = 0.0f;
    };

    sampleDirectLightFromFirstIntersection(
            Gbuffer, args.position, model, spp_direct,
            renderData, radiance_Dd, radiance_Ds);
    calcVar(radiance_Dd, exposure);
    calcVar(radiance_Ds, exposure);

    std::vector<LightSample> samples;
    for (int T = 0; T < spp_indirect; T++)
        sampleIndirectLightFromFirstIntersection(
                Gbuffer, args.position, base_diff, model,
                renderData, samples);
    if (samples.empty())
        return;

    mulWeight(samples, 1.0f / static_cast<float>(spp_indirect));

    float meanClum = 0.0f;
    for (const auto &sample: samples)
        meanClum += dot(sample.bsdfPdf * sample.light, RGB_Weight) * sample.weight;

    static const float clampThreshold = 16.0f;
    for (auto &sample: samples) {
        float Clum = dot(sample.bsdfPdf * sample.light, RGB_Weight) * sample.weight;
        if (Clum / (meanClum - Clum + eps_zero) > clampThreshold)
            continue;
        accumulateInwardRadiance(
                Gbuffer.opacity < eps_zero ? vec3(0.0f) : Gbuffer.baseColor,
                sample, radiance_Id, radiance_Is
        );
    }
    calcVar(radiance_Id, exposure);
    calcVar(radiance_Is, exposure);
}

struct PixelPos {
    int x, y;

    PixelPos(int x, int y) : x(x), y(y) {}
};

struct TaskQueue {
    std::vector<std::vector<PixelPos>> tasks;
    std::mutex mtx;

    int getTask(std::vector<PixelPos> &task) {
        std::lock_guard<std::mutex> lock(mtx);
        if (tasks.empty())
            return -1;
        task = tasks.back();
        tasks.pop_back();
        std::cout << "Task: " << tasks.size() << std::endl;
        return 0;
    }
};

void render(const Model &model, const RenderArgs &args,
            RenderData &renderData, Image &image, const std::vector<PixelPos> &pixels) {
    for (const auto &pixel: pixels)
        renderPixel(model, args, renderData, image, pixel.x, pixel.y);
}

void renderWorker(const Model &model, const RenderArgs &args,
                  TaskQueue &taskQueue, RenderData &renderData, Image &image) {
    while (true) {
        std::vector<PixelPos> pixels;
        int flag = taskQueue.getTask(pixels);
        if (flag == -1)
            break;
        render(model, args, renderData, image, pixels);
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
    TaskQueue taskQueue;
    for (int y = 0; y < height; y++) {
        std::vector<PixelPos> task;
        task.reserve(width);
        for (int x = 0; x < width; x++)
            task.emplace_back(x, y);
        taskQueue.tasks.push_back(task);
    }

    std::vector<RenderData> datas;
    datas.reserve(threads);
    for (int i = 0; i < threads; i++)
        datas.emplace_back(i);

    std::vector<std::thread> threadPool;
    threadPool.reserve(threads);
    for (int i = 0; i < threads; i++)
        threadPool.emplace_back([&](int index) {
            renderWorker(model, args, taskQueue, datas[index], image);
        }, i);
    for (auto &thread: threadPool)
        thread.join();

    int C_lightSamples = 0;
    for (const auto &data: datas)
        C_lightSamples += data.C_lightSamples;

    std::cout << "Rendering completed in " << clock() - startTime << " ms." << std::endl;
    std::cout << "Direct light samples: " << C_lightSamples << std::endl;

    auto exportImage = [&]
            (const std::string &tag, int shadeOptions, float exposure) {
        image.postProcessing(shadeOptions, exposure);
        image.save((args.savePath + "(" + tag + ").png").c_str());
    };

    exportImage("DiffuseColor", Image::BaseColor, exposure);
    exportImage("DiffuseColor_FXAA", Image::BaseColor | Image::DoFXAA, exposure);
    exportImage("shapeNormal", Image::shapeNormal, exposure);
    exportImage("surfaceNormal", Image::surfaceNormal, exposure);

    exportImage("Direct_Diffuse", Image::Direct_Diffuse, exposure);
    exportImage("Direct_Specular", Image::Direct_Specular, exposure);
    exportImage("Indirect_Diffuse", Image::Indirect_Diffuse, exposure);
    exportImage("Indirect_Specular", Image::Indirect_Specular, exposure);
    exportImage("Raw", Image::Full, exposure);
    exportImage("Raw_Bloom", Image::Full | Image::DoBloom, exposure);
    exportImage("Raw_FXAA", Image::Full | Image::DoFXAA, exposure);
    image.filter();
    exportImage("Direct_Diffuse_Filter", Image::Direct_Diffuse, exposure);
    exportImage("Direct_Specular_Filter", Image::Direct_Specular, exposure);
    exportImage("Indirect_Diffuse_Filter", Image::Indirect_Diffuse, exposure);
    exportImage("Indirect_Specular_Filter", Image::Indirect_Specular, exposure);
    exportImage("Filter", Image::Full, exposure);
    exportImage("Filter_Bloom", Image::Full | Image::DoBloom, exposure);
    exportImage("Filter_FXAA", Image::Full | Image::DoFXAA, exposure);

    std::cout << "Post processing finished. Total: " << clock() - startTime << " ms." << std::endl;
}