#include "sampling.h"
#include "disney/disney.h"

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

float GTR1(float NdotH, float a)
{
    if (a >= 1) return 1/PI;
    float a2 = a*a;
    float t = 1 + (a2-1)*NdotH*NdotH;
    return (a2-1) / (PI*log(a2)*t);
}

float GTR2(float NdotH, float a)
{
    float a2 = a*a;
    float t = 1 + (a2-1)*NdotH*NdotH;
    return a2 / (PI * t*t);
}

float GTR2_aniso(float NdotH, float HdotX, float HdotY, float ax, float ay)
{
    return 1 / (PI * ax*ay * sqr( sqr(HdotX/ax) + sqr(HdotY/ay) + NdotH*NdotH ));
}

float smithG_GGX(float NdotV, float alphaG)
{
    float a = alphaG*alphaG;
    float b = NdotV*NdotV;
    return 1 / (NdotV + sqrt(a + b - a*b));
}

float smithG_GGX_aniso(float NdotV, float VdotX, float VdotY, float ax, float ay)
{
    return 1 / (NdotV + sqrt( sqr(VdotX*ax) + sqr(VdotY*ay) + sqr(NdotV) ));
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
    if (NdotL < 0 || NdotV < 0) return vec3(0);

    vec3 H = normalize(L+V);
    float NdotH = dot(N,H);
    float LdotH = dot(L,H);

    vec3 Cdlin = surface.baseColor;
    float Cdlum = dot(Cdlin, RGB_Weight); // luminance approx.

    const float subsurface = 0.0f;
    const float specular = 0.0f;
    const float specularTint = 1.0f;
    const float sheen = 0.0f;
    const float sheenTint = 0.0f;
    const float anisotropic = 0.0f;
    const float clearcoat = 1.0f;
    const float clearcoatGloss = 0.2f;
    const float clearcoatTint = 0.0f;

    vec3 Ctint = Cdlum > 0 ? Cdlin/Cdlum : vec3(1); // normalize lum. to isolate hue+sat
    vec3 Cspec0 = mix(specular * .08f * glm::mix(vec3(1), Ctint, specularTint), Cdlin, surface.metallic);
    vec3 Csheen = mix(vec3(1), Ctint, sheenTint);

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
    float aspect = sqrt(1-anisotropic*.9);
    float ax = std::max(.001f, sqr(surface.roughness)/aspect);
    float ay = std::max(.001f, sqr(surface.roughness)*aspect);
    float Ds = GTR2_aniso(NdotH, dot(H, X), dot(H, Y), ax, ay);
    float FH = SchlickFresnel(LdotH);
    vec3 Fs = mix(Cspec0, vec3(1.0f), FH);
    float Gs;
    Gs  = smithG_GGX_aniso(NdotL, dot(L, X), dot(L, Y), ax, ay);
    Gs *= smithG_GGX_aniso(NdotV, dot(V, X), dot(V, Y), ax, ay);

    // sheen
    vec3 Fsheen = FH * sheen * Csheen;

    // clearcoat (ior = 1.5 -> F0 = 0.04)
    float Dr = GTR1(NdotH, mix(.1,.001,clearcoatGloss));
    float Fr = mix(.04, 1.0, FH);
    float Gr = smithG_GGX(NdotL, .25) * smithG_GGX(NdotV, .25);
    vec3 Vclearcoat = (0.25f*clearcoat*Gr*Fr*Dr) * (1 - surface.metallic)
            * glm::mix(vec3(1), Ctint, clearcoatTint);

    vec3 ret = ((1.0f/M_PI) * mix(Fd, ss, subsurface) * Cdlin + Fsheen) * (1-surface.metallic)
               + Vclearcoat
               + Gs*Fs*Ds;

    return ret;
}

vec3 BRDF::getBRDFWithCos(vec3 outDir) const {
    return getBRDF(outDir) * dot(outDir, surface.surfaceNormal);
}

void BRDF::sample(Generator &gen, vec3 &Dir, vec3 &brdfPdf) const {
    auto sample0 = [&]() -> vec3 {
        float u = gen(), v = gen() * 2.0f * M_PI;
        float d = sqrt(u);
        float z = sqrt(1 - d*d);
        float x = d * cos(v), y = d * sin(v);
        return x * tangent + y * bitangent + z * surface.surfaceNormal;
    };
    for (int T = 1; T <= MaxTrys; T++) {
        Dir = sample0();
        if (dot(Dir, surface.shapeNormal) > 0.0f) {
            brdfPdf = getBRDF(Dir);
            brdfPdf *= 1.0f / T;
        }
    }
    Dir = vec3(0.0f);
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

void sampleLightFace(const vec3 &pos, const LightObject &lightObject, Generator &gen, float &P_success, int &faceIndex) {
    vec3 lightDir;
    for(int T = 1; T <= MaxTrys; T++) {
        faceIndex = lightObject.faceDist(gen);
        auto &lightFace = lightObject.lightFaces[faceIndex];
        lightDir = normalize(lightFace.position - pos);
        float cosPhi = dot(lightFace.normal, -lightDir);
        if (cosPhi <= 0.0f || gen() > cosPhi) {
            if (T == MaxTrys) {
                faceIndex = -1;
                return;
            }
        }
        else {
            P_success = 1.0f / T;
            return;
        }
    }
}

vec3 sampleDirectLight(const vec3 &pos, const BRDF &brdf, const Model &model,
                       Generator &gen) {

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

    float P_success;
    int faceIndex;
    sampleLightFace(pos, lightObject, gen, P_success, faceIndex);
    if (faceIndex == -1)
        return vec3(0.0f);
    auto &lightFace = lightObject.lightFaces[faceIndex];

    vec3 lightDir = normalize(lightFace.position - pos);
    float distance = length(lightFace.position - pos);
    if (model.rayHit_test({pos, lightDir}, distance - eps_lightRadius))
        return vec3(0.0f);

    float distanceSquared = distance * distance + eps_lightRadius;
    float contribution = lightObject.power / (distanceSquared * P_lightObject) * P_success;
    return lightObject.color * contribution * brdf.getBRDF(lightDir);
}
