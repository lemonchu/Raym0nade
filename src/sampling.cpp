#include <iostream>
#include "sampling.h"

BRDF::BRDF(const vec3 &inDir) : inDir(inDir) {}

const float PI = 3.14159265358979323846;
static const vec3 RGB_Weight = vec3(0.3f, 0.6f, 0.1f);

float sqr(float x) { return x*x; }

float clamp(float x, float a, float b) {
    return x < a ? a : (x > b ? b : x);
}

float mix(float a, float b, float t) {
    return a*(1-t) + b*t;
}

float SchlickFresnel(float u) {
    float m = clamp(1-u, 0, 1);
    float m2 = m*m;
    return m2*m2*m; // pow(m,5)
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


void BRDF::genTangentSpace() {
    getTangentSpaceWithInDir(surface.surfaceNormal, inDir, tangent, bitangent);
}

vec3 BRDF::getBRDF(vec3 outDir) const
{
    const vec3 &V = inDir;
    const vec3 &L = outDir;
    const vec3 &N = surface.surfaceNormal;
    const vec3 &X = tangent;
    const vec3 &Y = bitangent;

    float NdotL = dot(N,L);
    float NdotV = dot(N,V);
    if (NdotL <= 0.0f || NdotV <= 0.0f)
        return vec3(0);

    vec3 H = normalize(L+V);
    float NdotH = dot(N,H);
    float LdotH = dot(L,H);

    vec3 Cdlin = surface.baseColor;
    float Cdlum = dot(Cdlin, RGB_Weight); // luminance approx.

    const float subsurface = 0.0f;
    const float specular = 0.1f;
    const float specularTint = 0.0f;
    const float sheen = 0.0f;
    const float sheenTint = 0.0f;
    const float clearcoat = 1.5f;
    const float clearcoatGloss = 0.25f;
    const float clearcoatTint = 0.0f;

    vec3 Ctint = Cdlum > 0.0f ? Cdlin/Cdlum : vec3(1); // normalize lum. to isolate hue+sat
    vec3 Cspec0 = mix(specular * .08f * glm::mix(vec3(1), Ctint, specularTint), Cdlin, surface.metallic);
    vec3 Csheen = mix(vec3(1.0f), Ctint, sheenTint);

    // Diffuse fresnel - go from 1 at normal incidence to .5 at grazing
    // and mix in diffuse retro-reflection based on roughness
    float FL = SchlickFresnel(NdotL), FV = SchlickFresnel(NdotV);
    float Fd90 = 0.5f + 2 * LdotH*LdotH * surface.roughness;
    float Fd = mix(1.0f, Fd90, FL) * mix(1.0f, Fd90, FV);

    // Based on Hanrahan-Krueger brdf approximation of isotropic bssrdf
    // 1.25 scale is used to (roughly) preserve albedo
    // Fss90 used to "flatten" retroreflection based on roughness
    float Fss90 = LdotH * LdotH * surface.roughness;
    float Fss = mix(1.0f, Fss90, FL) * mix(1.0f, Fss90, FV);
    float ss = 1.25f * (Fss * (1 / (NdotL + NdotV) - .5f) + .5f);

    // specular
    float Ds = GTR2(NdotH, surface.roughness);
    float FH = SchlickFresnel(LdotH);
    vec3 Fs = mix(Cspec0, vec3(1.0f), FH);
    float Gs;
    Gs  = smithG_GGX(NdotL, surface.roughness);
    Gs *= smithG_GGX(NdotV, surface.roughness);

    // sheen
    vec3 Fsheen = FH * sheen * Csheen;

    // clearcoat (ior = 1.5 -> F0 = 0.04)
    float Dr = GTR1(NdotH, mix(.1,.001,clearcoatGloss));
    float Fr = mix(.04, 1.0, FH);
    float Gr = smithG_GGX(NdotL, .25) * smithG_GGX(NdotV, .25);

    vec3 ret = (static_cast<float>(1.0f/M_PI) * mix(Fd, ss, subsurface) * Cdlin + Fsheen) * (1-surface.metallic)
             + (0.25f*clearcoat*Gr*Fr*Dr) * glm::mix(vec3(1), Ctint, clearcoatTint)
             + Gs*Fs*Ds;

    return ret * dot(outDir, surface.surfaceNormal);
}

static const int MaxTrys = 16;

void BRDF::sampleBRDF(Generator &gen, vec3 &outDir, vec3 &brdfPdf, int &fails) const {

    auto sampleGGX2 = [&](vec3 &outDir, float &pdf) -> void {
        float u = gen(), phi = gen() * 2.0f * M_PI;
        float
            cosTheta = sqrt((1.0f-u) / (1.0f+(sqr(surface.roughness)-1.0f)*u)),
            sinTheta = sqrt(1 - cosTheta*cosTheta);
        vec3 H = sinTheta * cos(phi) * tangent
                + sinTheta * sin(phi) * bitangent
                + cosTheta * surface.surfaceNormal;
        outDir = reflect(-inDir, H);
        float LDotH = dot(outDir, H);
        float LdotN = dot(outDir, surface.surfaceNormal);
        if (LDotH < 0.0f || LdotN < 0.0f) {
            pdf = 0.0f;
            return ;
        }
        pdf = GTR2(cosTheta, surface.roughness) * cosTheta / (4.0f * LDotH);
    };
    for (int T = 1; T <= MaxTrys; T++) {
        float pdf;
        sampleGGX2(outDir, pdf);
        if (dot(outDir, surface.shapeNormal) > 0.0f && pdf > 0.0f) {
            brdfPdf = getBRDF(outDir) / pdf;
            return;
        }
        fails++;
    }
    outDir = vec3(NAN);
}

void BRDF::sampleCos(Generator &gen, vec3 &outDir, vec3 &brdfPdf, int &fails) const {
    auto sampleCos = [&]() -> vec3 {
        float u = gen(), phi = gen() * 2.0f * M_PI;
        float d = sqrt(u);
        float z = sqrt(1 - d*d);
        float x = d * cos(phi), y = d * sin(phi);
        return x * tangent + y * bitangent + z * surface.surfaceNormal;
    };
    for (int T = 1; T <= MaxTrys; T++) {
        outDir = sampleCos();
        float pdf = dot(outDir, surface.surfaceNormal) / M_PI;
        if (dot(outDir, surface.shapeNormal) > 0.0f) {
            brdfPdf = getBRDF(outDir) / pdf;
            return;
        }
        fails++;
    }
    outDir = vec3(NAN);
}

void BRDF::sample(Generator &gen, vec3 &Dir, vec3 &brdfPdf, int &fails) const {
    if (gen() < surface.metallic)
        sampleBRDF(gen, Dir, brdfPdf, fails);
    else
        sampleCos(gen, Dir, brdfPdf, fails);
}

int sample(const std::vector<float> &weights, float randomValue) {
    float cumulativeWeight = 0.0f;
    for (int i = 0; i < weights.size(); ++i) {
        cumulativeWeight += weights[i];
        if (randomValue <= cumulativeWeight)
            return i;
    }
}

void sampleLightObject(const vec3 &pos, const BRDF &brdf, const Model &model,
                       Generator &gen, float &prob, int &lightIndex) {
    float totalWeight = 0.0f;
    std::vector<float> weights;
    for (const auto & lightObject : model.lightObjects)  {
        vec3 lightDir = normalize(lightObject.center - pos);
        float distance = glm::length(lightObject.center - pos);
        vec3 brdfPdf = brdf.getBRDF(lightDir);
        float albedo = dot((vec3)brdf.surface.baseColor, RGB_Weight);
        float weight = albedo * lightObject.power / (distance * distance + eps_lightRadius);
        weights.push_back(weight);
        totalWeight += weight;
    }
    if (totalWeight == 0.0f) {
        lightIndex = -1;
        return;
    }
    lightIndex = sample(weights, totalWeight * gen());
    prob = weights[lightIndex] / totalWeight;
}

void generateRandomPointInLightFace(const LightFace &lightFace, const vec3 &pos,
                                    Generator &gen, vec3 &lightPos, float &cosPhi) {
    float a = gen(), b = gen();
    if (a + b > 1.0f) {
        a = 1.0f - a;
        b = 1.0f - b;
    }
    float c = 1.0f - a - b;
    const Face &face = lightFace.face;
    lightPos = a * face.v[0] + b * face.v[1] + c * face.v[2];
    vec3 lightDir = normalize(lightPos - pos);
    vec2 texUV =
            a * face.data[0]->uv +
            b * face.data[1]->uv +
            c * face.data[2]->uv;
    vec3 shapeNormal, surfaceNormal;
    getHitNormal(face, lightDir, vec3(a,b,c), texUV, shapeNormal, surfaceNormal);
    cosPhi = dot(surfaceNormal, -lightDir);
}

void sampleLightFace(const vec3 &pos, const LightObject &lightObject,
                     Generator &gen, vec3 &lightPos, float &faceFactor, int &fails) {
    vec3 lightDir;
    static const int MaxTrys = 16;
    for (int T = 1; T <= MaxTrys; T++) {
        int faceIndex = lightObject.faceDist(gen);
        auto &lightFace = lightObject.lightFaces[faceIndex];
        float cosPhi;
        generateRandomPointInLightFace(lightFace, pos, gen, lightPos, cosPhi);
        if (cosPhi > 0.0f) {
            faceFactor = cosPhi / (dot(pos - lightPos, pos - lightPos) + eps_lightRadius);
            return;
        }
    }
    lightPos = vec3(NAN);
}

vec3 sampleDirectLight(const vec3 &pos, const BRDF &brdf, const Model &model,
                       Generator &gen, vec3 &brdfPdf, int &fails) {

    // vec3 light = vec3(0.0f);
    // auto &lightObjects = model.lightObjects;
    // for (const auto& lightObject : lightObjects) {
    //     float distance = glm::length(lightObject.center - pos);
    //     float dot = std::max(glm::dot(normalize(lightObject.center - pos), info.normal), 0.00f);
    //     light += dot * lightObject.color * lightObject.power / (distance * distance + eps_lightRadius);
    // }
    // return vec4(light, 1.0f);

    float P_lightObject;
    int lightIndex;
    sampleLightObject(pos, brdf,  model, gen, P_lightObject, lightIndex);
    if (lightIndex == -1)
        return vec3(0.0f);
    auto& lightObject = model.lightObjects[lightIndex];

    vec3 lightPos;
    float faceFactor;
    sampleLightFace(pos, lightObject, gen, lightPos, faceFactor, fails);
    if (isnan(lightPos))
        return vec3(0.0f);

    vec3 lightDir = normalize(lightPos - pos);
    float distance = length(lightPos - pos);
    if (model.rayHit_test({pos, lightDir}, distance - eps_lightRadius))
        return vec3(0.0f);

    float contribution = lightObject.power * faceFactor / P_lightObject;
    brdfPdf = brdf.getBRDF(lightDir);
    if (isnan(brdfPdf))
        std::cerr << "Warning: NaN detected! (sampleDirectLight::brdfPdf)" << std::endl;
    if (isnan(contribution))
        std::cerr << "Warning: NaN detected! (sampleDirectLight::contribution)" << std::endl;
    return lightObject.color * contribution;
}