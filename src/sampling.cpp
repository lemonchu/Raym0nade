#include <iostream>
#include "sampling.h"

BSDF::BSDF(const vec3 &inDir, const vec3 &hit_position) : inDir(inDir) {
    surface.position = hit_position;
}

BSDF::BSDF(const vec3 &inDir, const HitInfo &hitInfo) : inDir(inDir), surface(hitInfo) {}

float sqr(float x) { return x*x; }

float clamp(float x, float a, float b) {
    return x < a ? a : (x > b ? b : x);
}

float mix(float a, float b, float t) {
    return a*(1-t) + b*t;
}

float pow5(float x) {
    float x2 = x*x;
    return x2*x2*x;
}

float SchlickFresnel(float u) {
    float m = clamp(1-u, 0, 1);
    return pow5(m);
}

float GTR1(float NdotH, float a) {
    if (a >= 1) return 1/PI;
    float a2 = a*a;
    float t = 1 + (a2-1)*NdotH*NdotH;
    return (a2-1) / (PI*log(a2)*t);
}

float GTR2(float NdotH, float a) {
    float a2 = a*a;
    float t = 1 + (a2-1)*NdotH*NdotH;
    return a2 / (PI * t*t);
}

float smithG_GGX(float NdotV, float alphaG) {
    float a = alphaG*alphaG;
    float b = NdotV*NdotV;
    return 1 / (NdotV + sqrt(a + b - a*b));
}

vec3 BSDF::getBRDF(vec3 L) const {
    const vec3 &V = inDir;
    const vec3 &N = surface.surfaceNormal;

    float NdotL = dot(N,L);
    float NdotV = dot(N,V);
    if (NdotL <= 0.0f || NdotV <= 0.0f)
        return vec3(0);

    vec3 H = normalize(L+V);
    float NdotH = dot(N,H);
    float LdotH = dot(L,H);

    vec3 Cdlin = surface.baseColor;
    float Cdlum = dot(Cdlin, RGB_Weight); // luminance approx.

    const float
        subsurface = 0.0f,
        specularTint = 0.0f,
        sheen = 0.0f,
        sheenTint = 0.0f,
        clearcoat = 1.5f,
        clearcoatGloss = 0.2f,
        clearcoatTint = 0.0f;

    vec3 Ctint = Cdlum > 0.0f ? Cdlin/Cdlum : vec3(1.0f); // normalize lum. to isolate hue+sat
    vec3 Cspec0 = mix(surface.specular * glm::mix(vec3(1.0f), Ctint, specularTint), Cdlin, surface.metallic);
    vec3 Csheen = mix(vec3(1.0f), Ctint, sheenTint);

    // Diffuse fresnel - go from 1 at normal incidence to .5 at grazing
    // and mix in diffuse retro-reflection based on roughness
    float FL = SchlickFresnel(NdotL), FV = SchlickFresnel(NdotV);
    float Fd90 = 0.5f + 2.0f * LdotH*LdotH * surface.roughness;
    float Fd = mix(1.0f, Fd90, FL) * mix(1.0f, Fd90, FV);

    // Based on Hanrahan-Krueger brdf approximation of isotropic bssrdf
    // 1.25 scale is used to (roughly) preserve albedo
    // Fss90 used to "flatten" retroreflection based on roughness
    float Fss90 = LdotH * LdotH * surface.roughness;
    float Fss = mix(1.0f, Fss90, FL) * mix(1.0f, Fss90, FV);
    float ss = 1.25f * (Fss * (1.0f / (NdotL + NdotV) - .5f) + .5f);

    // specular
    float Ds = GTR2(NdotH, surface.roughness);
    float FH;
    if (surface.eta <= 1.0f) {
        FH = SchlickFresnel(LdotH);
    } else {
        float cosThetaI = fabs(LdotH);
        float sinThetaI = sqrt_s(1.0f - cosThetaI * cosThetaI);
        float sinThetaT = surface.eta * sinThetaI;
        if (sinThetaT >= 1.0f) {
            // Total internal reflection
            FH = 1.0f;
        } else {
            float cosThetaT = sqrt_s(1.0f - sinThetaT * sinThetaT);
            float R0 = (surface.eta - 1.0f) / (surface.eta + 1.0f);
            R0 = R0 * R0;
            FH = mix(R0, 1.0f, SchlickFresnel(cosThetaT));
        }
    }
    vec3 Fs = mix(Cspec0, vec3(1.0f), FH);
    float Gs;
    Gs  = smithG_GGX(NdotL, surface.roughness);
    Gs *= smithG_GGX(NdotV, surface.roughness);

    // sheen
    vec3 Fsheen = FH * sheen * Csheen;

    // clearcoat (ior = 1.5 -> F0 = 0.04)
    float Dr = GTR1(NdotH, mix(.1f,.001f,clearcoatGloss));
    float Fr = mix(.04f, 1.0f, FH);
    float Gr = smithG_GGX(NdotL, .25f) * smithG_GGX(NdotV, .25f);

    vec3 ret = ((1.0f/PI) * mix(Fd, ss, subsurface) * Cdlin + Fsheen);
    ret *= (1.0f - surface.metallic);
    ret += (0.25f*clearcoat*Gr*Fr*Dr) * glm::mix(vec3(1.0f), Ctint, clearcoatTint);
    ret *= surface.opacity;
    ret += Fs * Ds * Gs;

#ifdef RAY_DEBUG
    if (!isfinite(ret)) {
        std::cerr << "Wrong BRDF (getBRDF)" << std::endl;
        std::cerr << "ret: " << ret.x << " " << ret.y << " " << ret.z << std::endl;
        std::cerr << "Cdlin: " << Cdlin.x << " " << Cdlin.y << " " << Cdlin.z << std::endl;
        std::cerr << "Ctint: " << Ctint.x << " " << Ctint.y << " " << Ctint.z << std::endl;
        std::cerr << "Cspec0: " << Cspec0.x << " " << Cspec0.y << " " << Cspec0.z << std::endl;
        std::cerr << "Csheen: " << Csheen.x << " " << Csheen.y << " " << Csheen.z << std::endl;
        std::cerr << "Fs: " << Fs.x << " " << Fs.y << " " << Fs.z << std::endl;
        std::cerr << "Fsheen: " << Fsheen.x << " " << Fsheen.y << " " << Fsheen.z << std::endl;
        std::cerr << "Fd: " << Fd << std::endl;
        std::cerr << "ss: " << ss << std::endl;
        std::cerr << "Ds: " << Ds << std::endl;
        std::cerr << "FH: " << FH << std::endl;
        std::cerr << "Fs: " << Fs.x << " " << Fs.y << " " << Fs.z << std::endl;
        std::cerr << "Gs: " << Gs << std::endl;
        std::cerr << "Dr: " << Dr << std::endl;
        std::cerr << "Fr: " << Fr << std::endl;
        std::cerr << "Gr: " << Gr << std::endl;
        std::cerr << "surface.metallic: " << surface.metallic << std::endl;
        std::cerr << "surface.opacity: " << surface.opacity << std::endl;
        std::cerr << "surface.roughness: " << surface.roughness << std::endl;
        std::cerr << "surface.baseColor: " << surface.baseColor.x << " " << surface.baseColor.y << " " << surface.baseColor.z << std::endl;
        std::cerr << "surface.surfaceNormal: " << surface.surfaceNormal.x << " " << surface.surfaceNormal.y << " " << surface.surfaceNormal.z << std::endl;
    }
#endif

    return ret * NdotL;
}

vec3 BSDF::getBTDF(vec3 L) const {
    const vec3 &N = surface.surfaceNormal;
    float NdotL = dot(surface.surfaceNormal, L);
    if (NdotL >= 0.0f)
        return vec3(0.0f);

    const vec3 &V = inDir;

    vec3 H = normalize(L + surface.eta * V);

    float D = GTR2(dot(N, H), surface.roughness);

    float BTDF =  D  * (-NdotL);

    float LdotH = dot(L, H);
    float NdotV = dot(N, V);
    float HdotV = dot(H, V);
    BTDF *= abs(LdotH * HdotV)/(abs(NdotL * NdotV) + eps_zero);
    float k = surface.eta * HdotV + LdotH;
    BTDF /= k*k;

    vec3 ret = vec3(BTDF);

#ifdef RAY_DEBUG
    if (!isfinite(ret)) {
        std::cerr << "Wrong BTDF (getBTDF)" << std::endl;
        std::cerr << "BTDF: " << BTDF << std::endl;
        std::cerr << "k: " << k << std::endl;
        std::cerr << "D: " << D << std::endl;
    }
#endif

    return ret;
}

vec3 BSDF::getBSDF(vec3 outDir) const {
    if (surface.opacity > 1.0f - eps_zero || surface.entering) // 在非透明面上，或进入透明面时，不考虑直接光照
        return getBRDF(outDir);
    else
        return getBTDF(outDir);
}

static const int MaxTrys = 16;

void BSDF::sampleBRDF(Generator &gen, vec3 &outDir, vec3 &brdfPdf, float &pdf, int &fails) const {

    vec3 tangent, bitangent;
    getTangentSpaceWithInDir(surface.surfaceNormal, inDir, tangent, bitangent);

    auto sampleGTR2 = [&](vec3 &outDir, float &pdf) -> void {
        float u = gen(), phi = gen() * 2.0f * PI;
        float
            cosTheta = sqrt((1.0f-u) / (1.0f+(sqr(surface.roughness)-1.0f)*u)),
            sinTheta = sqrt(1 - cosTheta*cosTheta);
        vec3 H = sinTheta * cos(phi) * tangent
                + sinTheta * sin(phi) * bitangent
                + cosTheta * surface.surfaceNormal;
        outDir = 2.0f * dot(inDir, H) * H - inDir;
        float LDotH = dot(outDir, H);
        float LdotN = dot(outDir, surface.surfaceNormal);
        if (LDotH <= 0.0f || LdotN <= 0.0f) {
            pdf = 0.0f;
            return ;
        }
        pdf = GTR2(cosTheta, surface.roughness) / (4.0f * LDotH);
    };
    for (int T = 1; T <= MaxTrys; T++) {
        sampleGTR2(outDir, pdf);
        if (dot(outDir, surface.shapeNormal) > 0.0f && pdf > 0.0f) {
            brdfPdf = getBRDF(outDir) / pdf;
            return;
        }
        fails++;
    }
    pdf = 0.0f;
    brdfPdf = vec3(0.0f);
    outDir = vec3(NAN);
}

const float clampThreshold = 64.0f;
void clamp(vec3 &pdf) {
    float Clum = dot(pdf, RGB_Weight);
    if (Clum > clampThreshold)
        pdf /= Clum / clampThreshold;
}

void BSDF::sampleCos(Generator &gen, vec3 &outDir, vec3 &brdfPdf, float &pdf, int &fails) const {
    vec3 tangent, bitangent;
    getTangentSpaceWithInDir(surface.surfaceNormal, inDir, tangent, bitangent);

    auto sampleCos = [&]() -> vec3 {
        float u = gen(), phi = gen() * 2.0f * PI;
        float d = sqrt(u);
        float z = sqrt(1 - d*d);
        float x = d * cos(phi), y = d * sin(phi);
        return x * tangent + y * bitangent + z * surface.surfaceNormal;
    };
    for (int T = 1; T <= MaxTrys; T++) {
        outDir = sampleCos();
        pdf = dot(outDir, surface.surfaceNormal) / PI;
        if (dot(outDir, surface.shapeNormal) > 0.0f && pdf > 0.0f) {
            brdfPdf = getBRDF(outDir) / pdf;
            clamp(brdfPdf);
            return ;
        }
        fails++;
    }
    pdf = 0.0f;
    brdfPdf = vec3(0.0f);
    outDir = vec3(NAN);
}

void BSDF::preciseRefraction(vec3 &outDir, float &F) const {
    const vec3 &V = inDir;
    const vec3 &N = surface.surfaceNormal;
    float eta = surface.eta;

    float NdotV = dot(N, V);
    if (NdotV < 0.0f) {
        //std::cout << "        NdotV < 0.0f" << std::endl;
        F = 1.0f;
        outDir = vec3(NAN);
        return;
    }
    float delta = 1.0f - eta*eta*(1.0f - NdotV*NdotV);
    if (delta < 0.0f) {
        //std::cout << "        delta < 0.0f" << std::endl;
        F = 1.0f;
        outDir = vec3(NAN);
        return;
    }

    float k = eta * NdotV - sqrt(delta);
    outDir = k * N - eta * V;

    float LdotN = -dot(outDir, N);
    // Fresnel term using Schlick's approximation
    float R0 = (1.0f - eta) / (1.0f + eta);
    R0 = R0 * R0;
    F = R0 + (1.0f - R0) * SchlickFresnel(LdotN);

    /*if (rand()%4096 == 0) {
        std::cout << "surface.eta: " << surface.eta << std::endl;
        std::cout << "outDir: " << outDir.x << " " << outDir.y << " " << outDir.z << std::endl;
        std::cout << "V: " << V.x << " " << V.y << " " << V.z << std::endl;
        std::cout << "N: " << N.x << " " << N.y << " " << N.z << std::endl;
        std::cout << "dot(N, outDir): " << dot(N, outDir) << std::endl;
        std::cout << "dot(N, V): " << dot(N, V) << std::endl;
    }*/
}

void BSDF::sampleReflection(Generator &gen, vec3 &Dir, vec3 &brdfPdf, int &fails) const {

    vec3 Dir1, Dir2, brdfPdf1, brdfPdf2;
    int fail1 = 0, fail2 = 0;
    float pdf1, pdf2;
    sampleCos(gen, Dir1, brdfPdf1, pdf1, fail1);
    sampleBRDF(gen, Dir2, brdfPdf2, pdf2, fail2);
    float p1 = pdf1 / (pdf1 + pdf2);
/*
    std::cout << "             roughness: " << surface.roughness << std::endl;
    std::cout << "             BaseColor: " << surface.baseColor.x << " " << surface.baseColor.y << " " << surface.baseColor.z << std::endl;
    std::cout << "             surfaceNormal: " << surface.surfaceNormal.x << " " << surface.surfaceNormal.y << " " << surface.surfaceNormal.z << std::endl;
    std::cout << "             inDir: " << inDir.x << " " << inDir.y << " " << inDir.z << std::endl;
    std::cout << "             VdotN: " << dot(surface.surfaceNormal, inDir) << std::endl;
    std::cout << "             opacity: " << surface.opacity << std::endl;
    std::cout << "             eta: " << surface.eta << std::endl;
    std::cout << "             Dir1: " << Dir1.x << " " << Dir1.y << " " << Dir1.z << std::endl;
    std::cout << "             Dir2: " << Dir2.x << " " << Dir2.y << " " << Dir2.z << std::endl;
    std::cout << "             brdfPdf1: " << brdfPdf1.x << " " << brdfPdf1.y << " " << brdfPdf1.z << std::endl;
    std::cout << "             brdfPdf2: " << brdfPdf2.x << " " << brdfPdf2.y << " " << brdfPdf2.z << std::endl;
*/
    if (gen() < p1) { // 多重重要性采样 (Multiple Importance Sampling)
        Dir = Dir1;
        brdfPdf = brdfPdf1;
        fails += fail1;
    } else {
        Dir = Dir2;
        brdfPdf = brdfPdf2;
        fails += fail2;
    }
}

vec3 refract(const vec3 &V, const vec3 &N, float eta) {
    float NdotV = dot(N, V);
    float delta = 1.0f - eta * eta * (1.0f - NdotV * NdotV);
    if (delta < 0.0f)
        return vec3(NAN);
    float k = eta * NdotV - sqrt(delta);
    return k * N - eta * V;
}

void BSDF::sampleBTDF(Generator &gen, vec3 &outDir, vec3 &btdfPdf, int &fails) const {
    if (abs(surface.eta - 1.0f) < eps_zero) {
        outDir = -inDir;
        btdfPdf = vec3(1.0f);
        return;
    }
    vec3 tangent, bitangent;
    getTangentSpaceWithInDir(surface.surfaceNormal, inDir, tangent, bitangent);
    const vec3 &V = inDir;

    auto sampleGTR2 = [&](vec3 &outDir, float &weight) -> void {
        float u = gen(), phi = gen() * 2.0f * PI;
        float
                cosTheta = sqrt((1.0f-u) / (1.0f+(sqr(surface.roughness)-1.0f)*u)),
                sinTheta = sqrt(1 - cosTheta*cosTheta);
        vec3 H = sinTheta * cos(phi) * tangent
                 + sinTheta * sin(phi) * bitangent
                 + cosTheta * surface.surfaceNormal;
        outDir = refract(V, H, surface.eta);

        if (!isfinite(outDir)) {
            weight = 0.0f;
            return;
        }
        float VdotH = dot(V, H);
        float VdotN = dot(V, surface.surfaceNormal);
        float HdotN = dot(H, surface.surfaceNormal);
        weight = abs(VdotH) / abs(VdotN * HdotN);
    };
    for (int T = 1; T <= MaxTrys; T++) {
        float weight;
        sampleGTR2(outDir, weight);
        if (weight > 0.0f && dot(outDir, surface.shapeNormal) < 0.0f) {
            btdfPdf = vec3(weight);
            return;
        }
        fails++;
    }
    btdfPdf = vec3(0.0f);
    outDir = vec3(NAN);
}

LightSample::LightSample(const vec3 &bsdfPdf, const vec3 &light, float weight) :
        bsdfPdf(bsdfPdf), light(light), weight(weight) {}

int sample(const std::vector<float> &weights, float randomValue) {
    float cumulativeWeight = 0.0f;
    for (int i = 0; i < weights.size(); ++i) {
        cumulativeWeight += weights[i];
        if (randomValue <= cumulativeWeight)
            return i;
    }
    return -1;
}

void getLightObjectWeight(const vec3 &pos, const BSDF &bsdf, const Model &model,
                          std::vector<float> &weights) {
    weights.reserve(model.lightObjects.size());
    for (const auto & lightObject : model.lightObjects)  {
        vec3 lightDir = normalize(lightObject.center - pos);
        float distance = glm::length(lightObject.center - pos);
        vec3 bsdfPdf = bsdf.getBSDF(lightDir);
        float Clum = dot(bsdfPdf, RGB_Weight);
        float weight = Clum * lightObject.power / (distance * distance + eps_lightRadius);
        weights.emplace_back(weight);
    }
}

void generateRandomPointInLightFace(const Face &face, const vec3 &pos,
                                    Generator &gen, vec3 &lightPos, float &cosPhi) {
    float a = gen(), b = gen();
    if (a + b > 1.0f) {
        a = 1.0f - a;
        b = 1.0f - b;
    }
    float c = 1.0f - a - b;
    lightPos = a * face.v[0] + b * face.v[1] + c * face.v[2];
    vec3 lightDir = normalize(lightPos - pos);
    vec3 shapeNormal, surfaceNormal;
    bool entering;
    getHitNormals(face, lightDir, vec3(a, b, c),
                  shapeNormal, surfaceNormal, entering);
    cosPhi = dot(surfaceNormal, -lightDir);
}

void sampleLightFace(const vec3 &pos, const LightObject &lightObject,
                     Generator &gen, vec3 &lightPos, int &fails) {
    for (int T = 1; T <= MaxTrys; T++) {
        int faceIndex = lightObject.faceDist(gen);
        const Face& lightFace = lightObject.faces[faceIndex];
        float cosPhi;
        generateRandomPointInLightFace(lightFace, pos, gen, lightPos, cosPhi);
        if (cosPhi > 0.0f && gen() < cosPhi)
            return;
        fails++;
    }
    lightPos = vec3(NAN);
}

void sampleSkyBox(const SkyBox &skyBox, const vec3 &shapeNormal,
                  Generator &gen, vec3 &Dir, vec3 &light) {
    for (int T = 1; T <= MaxTrys; T++) {
        int pixelIndex = skyBox.dist(gen);
        int u = pixelIndex % skyBox.width, v = pixelIndex / skyBox.width;
        float phi = PI * (float(v)+0.5f) / float(skyBox.height);
        float theta = 2.0f * PI * (float(u)+0.5f) / float(skyBox.width);
        Dir = vec3(-sin(phi) * sin(theta), cos(phi), sin(phi) * cos(theta));
        float cosPhi = dot(Dir, shapeNormal);
        if (cosPhi > 0.0f) {
            light = skyBox.data[pixelIndex] / skyBox.dist.pdf(pixelIndex);
            return;
        }
    }
    Dir = vec3(NAN);
}

std::vector<LightSample> sampleDirectLight(
        const BSDF &bsdf, const Model &model, Generator &gen, int sampleCnt) {

    const vec3 &pos = bsdf.surface.position;
    std::vector<LightSample> samples;

    if (model.skyMap.empty()) {
        std::vector<float> weights;
        getLightObjectWeight(pos, bsdf, model, weights);
        float totalWeight = 0.0f;
        for (float weight : weights)
            totalWeight += weight;
        if (totalWeight == 0.0f)
            return {};

        for (int T = 0; T < sampleCnt; T++) {
            int lightIndex = sample(weights, totalWeight * gen());
            if (lightIndex == -1)
                continue;

            float P_lightObject = weights[lightIndex] / totalWeight;
            auto& lightObject = model.lightObjects[lightIndex];

            vec3 lightPos;
            int fails = 0;
            sampleLightFace(pos, lightObject, gen, lightPos, fails);
            if (!isfinite(lightPos))
                continue;

            vec3 lightDir = normalize(lightPos - pos);
            float distance = length(lightPos - pos);
            if (model.rayHit_test({pos, lightDir}, distance - eps_zero))
                continue;

            vec3 bsdfPdf = bsdf.getBSDF(lightDir);
            clamp(bsdfPdf);
            vec3 light = 2.0f * PI * lightObject.power * lightObject.color / (distance * distance + eps_lightRadius) / P_lightObject;

            if (!isfinite(bsdfPdf)) {
                std::cerr << "Wrong bsdf (sampleDirectLight_)" << std::endl;
                std::cerr << "lightDir: " << lightDir.x << " " << lightDir.y << " " << lightDir.z << std::endl;
                std::cerr << "bsdfPdf: " << bsdfPdf.x << " " << bsdfPdf.y << " " << bsdfPdf.z << std::endl;
                std::cerr << "P_lightObject: " << P_lightObject << std::endl;
            }
            samples.emplace_back(bsdfPdf, light, 1.0f / static_cast<float>(sampleCnt*(fails+1)));
        }
    } else {
        for (int T = 0; T < sampleCnt; T++) {
            vec3 Dir, light;
            sampleSkyBox(model.skyMap, bsdf.surface.surfaceNormal, gen, Dir, light);
            if (!isfinite(Dir))
                continue;
            if (model.rayHit_test({pos, Dir}, INFINITY))
                continue;
            vec3 bsdfPdf = bsdf.getBSDF(Dir);
            clamp(bsdfPdf);
            samples.emplace_back(bsdfPdf, light, 1.0f / static_cast<float>(sampleCnt));
        }
    }
    return samples;
}