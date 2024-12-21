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
    for (const auto &medium : mediums)
        ior = std::max(ior, medium.second.ior);
    return ior;
}

vec3 Medium::absorb() const {
    vec3 absorb = vec3(1.0f);
    for (const auto &medium : mediums)
        absorb *= medium.second.absorb;
    return absorb;
}

int Medium::size() const {
    return mediums.size();
}

void calc_dPdxy(const Ray &ray, float hit_t, const vec3 &normal, const RayDifferential &base_diff,
               vec3 &hit_dPdx, vec3 &hit_dPdy) {
    float dtdx = - dot(base_diff.dPdx + hit_t * base_diff.dDdx, normal) / dot(ray.direction, normal);
    float dtdy = - dot(base_diff.dPdy + hit_t * base_diff.dDdy, normal) / dot(ray.direction, normal);
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

    // Initialize the hit information Stage 2, we obtain the texture and the tangent space.
    calc_dPdxy(ray, hit.t_max, hitInfo.shapeNormal, base_diff, hit_dPdx, hit_dPdy);
    getHitMaterial(face, baryCoords, hit_dPdx, hit_dPdy, hitInfo);
}

// #define DEBUG_sampleRay

const float regularizationFactor = 0.24f;

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

static float P0_reflect = 0.24f;

vec3 sampleRay(Ray ray, const RayDifferential &base_diff, const Model &model,
               RenderData &renderData, int &fails, float roughnessFactor, int depth) {

    static const int maxRayDepth = 16;

    HitRecord hit = model.rayHit(ray);
    if (hit.t_max == INFINITY)
        return vec3(0.0f); // No intersection

    // 获取表面信息
    BSDF bsdf(-ray.direction, ray.origin + ray.direction * hit.t_max);
    vec3 hit_dPdx, hit_dPdy;
    getHitInfo(hit, ray, base_diff, hit_dPdx, hit_dPdy, bsdf.surface);

    float ior = bsdf.surface.eta; // 保存绝对折射率

    // 路径正则化
    roughnessFactor = std::max(roughnessFactor, regularizationFactor * bsdf.surface.roughness);
    bsdf.surface.roughness = std::max(bsdf.surface.roughness, roughnessFactor);

    // 根据光滑程度决定 RR 的截止概率
    static const float P_RR_Rough = 0.5f, P_RR_Smooth = 0.85f;
    const float P_RR = P_RR_Smooth + (P_RR_Rough - P_RR_Smooth) * sqrt(bsdf.surface.roughness);

    // 不计直接光照
    if (length(bsdf.surface.emission) > eps_zero)
        return vec3(0.0f);

#ifdef DEBUG_sampleRay
    std::cout << "depth: " << depth << std::endl;
    std::cout << "entering: " << (bsdf.surface.entering ? "true" : "false") << std::endl ;
    std::cout << "name: " << hit.face->material->name << std::endl << std::endl;;
#endif

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
#ifdef DEBUG_sampleRay
        std::cout << "    eta: " << bsdf.surface.eta << std::endl;
        std::cout << "    refractDir: " << refractDir.x << " " << refractDir.y << " " << refractDir.z << std::endl;
        std::cout << "    P_reflect: " << P_reflect << std::endl;
#endif
    }

    // 计算介质的吸收率
    vec3 absorb = getAbsorb(renderData.mediums.absorb(), hit.t_max);
#ifdef DEBUG_sampleRay
    std::cout << "        Absorb_depth: " << hit.t_max << std::endl;
    std::cout << "        Absorb: " << absorb.x << " " << absorb.y << " " << absorb.z << std::endl;
#endif

    vec3 newDir, bsdfPdf(NAN);
    RayDifferential next_diff;

    if (renderData.gen() <= P_reflect) {
        // 发生反射
#ifdef DEBUG_sampleRay
        std::cout << "reflect" << std::endl;
#endif
        bool doDirectLightSample =
                bsdf.surface.entering &&
                renderData.mediums.size() == 1;
        // 仅在空气中进行直接光照采样（在空气中反射）

        if (doDirectLightSample) {
            if (depth >= maxRayDepth || renderData.gen() > P_RR) {
                renderData.C_lightSamples++;
                vec3 light = sampleDirectLight(bsdf, model, renderData.gen, bsdfPdf);
                bsdfPdf /= P_reflect;
                if (depth < maxRayDepth)
                    bsdfPdf /= (1.0f - P_RR);
                if (!isfinite(light)) {
                    std::cout << "light is NAN! (sampleDirectLight)" << std::endl;
                }
                if (!isfinite(bsdfPdf)) {
                    std::cout << "bsdfPdf is NaN! (sampleDirectLight)" << std::endl;
                }
#ifdef DEBUG_sampleRay
                std::cout << "Light : " << light.x << " " << light.y << " " << light.z << std::endl;
                std::cout << "bsdfPdf : " << bsdfPdf.x << " " << bsdfPdf.y << " " << bsdfPdf.z << std::endl;
                std::cout << "depth :" << depth << std::endl;
                std::cout << std::endl;
#endif
                return bsdfPdf * light;
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
        bool doDirectLightSample =
                !bsdf.surface.entering &&
                renderData.mediums.size() == 2;
        // 仅在空气中进行直接光照采样（折射离开介质，进入空气）

        if (doDirectLightSample) {
            if (depth >= maxRayDepth || renderData.gen() > P_RR) {
                renderData.C_lightSamples++;
                vec3 light = sampleDirectLight(bsdf, model, renderData.gen, bsdfPdf);
                bsdfPdf /= (1.0f - P_reflect);
                if (depth < maxRayDepth)
                    bsdfPdf /= (1.0f - P_RR);
                if (!isfinite(light)) {
                    std::cout << "light is NAN! (sampleDirectLight.refract)" << std::endl;
                }
                if (!isfinite(bsdfPdf)) {
                    std::cout << "btdfPdf is NaN! (sampleDirectLight.refract)" << std::endl;
                }
#ifdef DEBUG_sampleRay
                    std::cout << "Light : " << light.x << " " << light.y << " " << light.z << std::endl;
                    std::cout << "bsdfPdf : " << bsdfPdf.x << " " << bsdfPdf.y << " " << bsdfPdf.z << std::endl;
                    std::cout << "depth :" << depth << std::endl;
                    std::cout << std::endl;
#endif
                return absorb * bsdfPdf * light;
            }
        }

        // 取折射方向（完美折射）
        newDir = refractDir;
        bsdfPdf = vec3((1.0f - F) * bsdf.surface.eta * bsdf.surface.eta);
        bsdfPdf /= (1.0f - P_reflect);
        if (doDirectLightSample)
            bsdfPdf /= P_RR;

        // Todo: 折射的 RayDifferential 计算
        next_diff = base_diff;
        if (!isfinite(bsdfPdf)) {
            std::cout << "bsdfPdf is NaN! (refract)" << P_reflect <<std::endl;
        }

        if (bsdf.surface.entering)
            renderData.mediums.insert(bsdf.surface.id, ior, bsdf.surface.baseColor);
        else
            renderData.mediums.erase(bsdf.surface.id);
    }

    // 若方向无效或递归过深，返回 0
    if (!isfinite(newDir) || depth == maxRayDepth)
        return vec3(0.0f); // Invalid direction, then return black

    Ray newRay = {bsdf.surface.position, newDir};

    bsdfPdf *= absorb;

#ifdef DEBUG_sampleRay
    std::cout << "newPos: " << newRay.origin.x << " " << newRay.origin.y << " " << newRay.origin.z << std::endl;
    std::cout << "newDir: " << newDir.x << " " << newDir.y << " " << newDir.z << std::endl;
    std::cout << "bsdfPdf: " << bsdfPdf.x << " " << bsdfPdf.y << " " << bsdfPdf.z << std::endl;
#endif

    vec3 irradiance = sampleRay(newRay, next_diff, model,
                                renderData, fails,
                                roughnessFactor, depth+1);

    if (!isfinite(irradiance)) {
        std::cout << "irradiance is NAN!" << std::endl;
    }
    if (!isfinite(bsdfPdf)) {
        std::cout << "bsdfPdf is NaN! " << depth << std::endl;
    }
    return  bsdfPdf * irradiance;
}

struct LightSample {
    vec3 bsdfPdf, light;
};

void sampleIndirectLightFromFirstIntersection(
        const HitInfo &hitInfo, const vec3 &origin, const RayDifferential &base_diff, const Model &model,
        RenderData &renderData, std::vector<LightSample> &samples, int &fails) {

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
        newDir = refractDir;

        bsdfPdf = vec3((1.0f - F) *bsdf.surface.eta * bsdf.surface.eta);
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

    vec3 light = sampleRay(newRay, next_diff, model,
                           renderData, fails,
                           roughnessFactor, 1);

#ifdef DEBUG_sampleRay
    std::cout << "Light : " << light.x << " " << light.y << " " << light.z << std::endl;
    std::cout << std::endl << std::endl;
#endif

    if (!isfinite(light)) {
        std::cerr << "light is NAN! (final)" << std::endl;
    }

    if (!isfinite(bsdfPdf)) {
        std::cerr << "bsdfPdf is NAN! (final)" << std::endl;
    }

    samples.push_back({bsdfPdf, light});
}

void accumulateInwardRadiance(RadianceData &radiance, vec3 inradiance) {
    if (!isfinite(inradiance)) {
        std::cerr << "NaN detected!" << std::endl;
        throw std::runtime_error("NaN detected!");
    }
    radiance.radiance += inradiance;
    radiance.Var += dot(inradiance, inradiance);
}

void accumulateInwardRadiance(const vec3 baseColor, const vec3 brdfPdf, const vec3 &light,
                              RadianceData &radiance_d, RadianceData &radiance_s) {

    if (light == vec3(0.0f))
        return;
    if (baseColor == vec3(0.0f)) {
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

void sampleDirectLightFromFirstIntersection(const HitInfo &hitInfo, const vec3 &origin, const Model &model,
                                            RenderData &renderData, RadianceData &radiance_Dd, RadianceData &radiance_Ds) {

    if (!isfinite(hitInfo.position) || length(hitInfo.emission) > 0)
        return; // 暂不计入直接光照

    vec3 inDir = normalize(hitInfo.position - origin);
    BSDF bsdf(-inDir, hitInfo);
    vec3 bsdfPdf = vec3(NAN);
    vec3 light = sampleDirectLight(bsdf, model, renderData.gen, bsdfPdf);
    accumulateInwardRadiance(hitInfo.baseColor, bsdfPdf, light, radiance_Dd, radiance_Ds);
}

void initRayDiff(const vec3 &d, const RenderArgs &args, RayDifferential &base_diff) {
    base_diff.dPdx = vec3(0.0f), base_diff.dPdy = vec3(0.0f);
    const vec3 dddx = args.accuracy * args.right, dddy = args.accuracy * args.up;
    float d_dot_d = dot(d, d), d_dot_dddx = dot(d, dddx), d_dot_dddy = dot(d, dddy);

    base_diff.dDdx = (d_dot_d * dddx - d * d_dot_dddx) / (sqrt(d_dot_d) * d_dot_d);
    base_diff.dDdy = (d_dot_d * dddy - d * d_dot_dddy) / (sqrt(d_dot_d) * d_dot_d);
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
                 RenderData &renderData, Image &image, int x, int y) {
    //renderData.gen.mt.seed(x^y);
#ifdef DEBUG_sampleRay
    if (x!=48||y!=41)
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
            rayX = static_cast<float>(x) - width / 2.0f,
            rayY = static_cast<float>(y) - height / 2.0f;
    vec3 d = view + accuracy * (rayX * right + rayY * up);
    HitInfo &Gbuffer = image.Gbuffer[id];

    RayDifferential base_diff;
    initRayDiff(d, args, base_diff);

    vec3 Dir = normalize(d);
    Ray ray = {position, Dir};

    HitRecord hit = model.rayHit(ray);
    if (hit.t_max == INFINITY) {
        Gbuffer.position = vec3(NAN);
        return;
    }

    Gbuffer.position = ray.origin + hit.t_max * ray.direction;
    vec3 hit_dPdx, hit_dPdy;
    getHitInfo(hit, ray, base_diff,hit_dPdx, hit_dPdy, Gbuffer);
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

    for (int T = 0; T < spp_direct; T++)
        sampleDirectLightFromFirstIntersection(
                Gbuffer, args.position, model,
                renderData, radiance_Dd, radiance_Ds);
    normalize(radiance_Dd, exposure, spp_direct);
    normalize(radiance_Ds, exposure, spp_direct);

    std::vector<LightSample> samples;
    int trys_indirect = 0;
    for (int T = 0; T < spp_indirect; T++)
        sampleIndirectLightFromFirstIntersection(
                Gbuffer, args.position, base_diff, model,
                renderData, samples, trys_indirect);
    if (samples.empty())
        return ;
    float meanClum = 0.0f;
    for (const auto &sample : samples)
        meanClum += dot(sample.bsdfPdf*sample.light, RGB_Weight);
    meanClum /= samples.size();

    static const float clampThreshold = 16.0f;
    for (auto &sample : samples) {
        float Clum = dot(sample.bsdfPdf*sample.light, RGB_Weight);
        float restClum = meanClum - Clum / samples.size();
        if (Clum/restClum > clampThreshold * samples.size())
            continue;
        accumulateInwardRadiance(Gbuffer.baseColor, sample.bsdfPdf, sample.light, radiance_Id, radiance_Is);
    }
    trys_indirect += spp_indirect;
    normalize(radiance_Id, exposure, trys_indirect);
    normalize(radiance_Is, exposure, trys_indirect);
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

/*
 * Render Workers
 * Master call this worker function to render the scene.
 */

void render(const Model &model, const RenderArgs &args,
            RenderData &renderData, Image &image, const std::vector<PixelPos> &pixels) {
    for (const auto &pixel : pixels)
        renderPixel(model, args, renderData, image, pixel.x, pixel.y);
}

void renderWorker(const Model &model, const RenderArgs &args,
                  TaskQueue &taskQueue, RenderData &renderData, Image &image) {
    while(true) {
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
        for (int x = 0; x < width; x++)
            task.emplace_back(x, y);
        taskQueue.tasks.push_back(task);
    }

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
            renderWorker(model, args, taskQueue, datas[index], image);
        }, i);
    for (auto &thread : threadPool)
        thread.join();

    int C_lightSamples = 0;
    for (const auto &data : datas)
        C_lightSamples += data.C_lightSamples;

    std::cout << "Rendering completed in " << clock() - startTime << " ms." << std::endl;
    std::cout << "Direct light samples: " << C_lightSamples << std::endl;

    auto exportImage = [&](const std::string &tag, int shadeOptions, bool exposure = true) {
        image.shade(exposure, shadeOptions);
        image.FXAA();
        image.gammaCorrection();
        if (exposure)
            image.save((args.savePath + "(" + tag + ")_no_FXAA.png").c_str());
        else
            image.save((args.savePath + "(" + tag + ").png").c_str());
    };

    exportImage("DiffuseColor", Image::BaseColor);
    exportImage("shapeNormal", Image::shapeNormal);
    exportImage("surfaceNormal", Image::surfaceNormal);

    exportImage("Direct_Diffuse", Image::Direct_Diffuse);
    exportImage("Direct_Specular", Image::Direct_Specular);
    exportImage("Indirect_Diffuse", Image::Indirect_Diffuse);
    exportImage("Indirect_Specular", Image::Indirect_Specular);
    exportImage("Raw", Image::Full);
    exportImage("Raw", Image::Full, false);
    image.filter();
    exportImage("Direct_Diffuse_Filter", Image::Direct_Diffuse);
    exportImage("Direct_Specular_Filter", Image::Direct_Specular);
    exportImage("Indirect_Diffuse_Filter", Image::Indirect_Diffuse);
    exportImage("Indirect_Specular_Filter", Image::Indirect_Specular);
    exportImage("Filter", Image::Full);
    exportImage("Filter", Image::Full, false);

    std::cout << "Post processing finished. Total: " << clock() - startTime << " ms." << std::endl;
}