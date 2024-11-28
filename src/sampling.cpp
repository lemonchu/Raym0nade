#include "sampling.h"

const float M_PI = 3.14159265359f;

BRDF::BRDF(const vec3 &inDir, const vec3 &shapeNormal, const vec3 &surfaceNormal) :
        inDir(inDir), shapeNormal(shapeNormal), surfaceNormal(surfaceNormal), roughness(1.0f), max(1.0f) {
    getTangentSpace(surfaceNormal, tangent, bitangent);
}

float BRDF::pdf(const vec3 &outDir) const {
    return dot(outDir,surfaceNormal) <= 0.0f ? 0.0f : std::max(dot(outDir, surfaceNormal), 0.0f);
}

float BRDF::P_accept(const vec3 &outDir) const {
    return pdf(outDir) / max;
}

vec3 BRDF::sample(std::mt19937 &gen) const {
    auto sample0 = [&]() -> vec3 {
        std::uniform_real_distribution<float> U(0.0f, 1.0f);
        float u = U(gen), v = U(gen) * 2.0f * M_PI;
        float d = sqrt(u);
        float z = sqrt(1 - d*d);
        float x = d * cos(v), y = d * sin(v);
        return x * tangent + y * bitangent + z * surfaceNormal;
    };
    vec3 direction;
    for (int T = 0; T < maxTrys; T++) {
        direction = sample0();
        if (dot(direction, shapeNormal) > 0.0f)
            break;
    }
    return direction;
}

void sampleLightObject(const vec3 &pos, const BRDF &brdf, const Model &model, std::mt19937 &gen, float &prob, int &lightIndex) {
    float totalWeight = 0.0f;
    std::vector<float> weights;
    for (const auto& lightObject : model.lightObjects) {
        float distance = glm::length(lightObject.center - pos);
        float dot = std::max(brdf.P_accept(normalize(lightObject.center - pos)) + 0.1f, 0.0f);
        float weight = dot * lightObject.power / (distance * distance + eps_lightRadius);
        weights.push_back(weight);
        totalWeight += weight;
    }
    if (totalWeight == 0.0f) {
        lightIndex = -1;
        return ;
    }
    lightIndex = std::discrete_distribution<int>(weights.begin(), weights.end())(gen);
    prob = weights[lightIndex] / totalWeight;
}

void sampleLightFace(const vec3 &pos, const LightObject &lightObject, std::mt19937 &gen, int &failCount, int &faceIndex) {
    vec3 lightDir;
    for(int T = 0;; T++) {
        faceIndex = lightObject.faceDist(gen);
        auto &lightFace = lightObject.lightFaces[faceIndex];
        lightDir = normalize(lightFace.position - pos);
        float cosPhi = dot(lightFace.normal, -lightDir);
        if (cosPhi <= 0.0f || std::uniform_real_distribution<float>(0.0f, 1.0f)(gen) > cosPhi) {
            if (T == maxTrys) {
                faceIndex = -1;
                return;
            }
            failCount++;
        }
        else
            break;
    }
}

void sampleDirectLight(const vec3 &pos, const BRDF &brdf, const Model &model, std::mt19937 &gen, int &failCount, vec3 &light) {

    // vec3 light = vec3(0.0f);
    // auto &lightObjects = model.lightObjects;
    // for (const auto& lightObject : lightObjects) {
    //     float distance = glm::length(lightObject.center - pos);
    //     float dot = std::max(glm::dot(normalize(lightObject.center - pos), info.normal), 0.00f);
    //     light += dot * lightObject.color * lightObject.power / (distance * distance + eps_lightRadius);
    // }
    // return vec4(light, 1.0f);

    light = vec3(0.0f);

    float P_lightObject;
    int lightIndex;
    sampleLightObject(pos, brdf, model, gen, P_lightObject, lightIndex);
    if (lightIndex == -1) return;
    auto& lightObject = model.lightObjects[lightIndex];

    int faceIndex;
    sampleLightFace(pos, lightObject, gen, failCount, faceIndex);
    if (faceIndex == -1) return;
    auto &lightFace = lightObject.lightFaces[faceIndex];

    vec3 lightDir = normalize(lightFace.position - pos);
    float distance = glm::length(lightFace.position - pos);
    if (model.rayHit_test({pos, lightDir}, distance - eps_lightRadius))
        return;

    float distanceSquared = distance * distance + eps_lightRadius;
    float contribution = lightObject.power * brdf.pdf(lightDir) / (distanceSquared * P_lightObject);
    light = lightObject.color * contribution;
}

Probe::Probe() {
    sampleCount = 0;
    normal = vec3(0);
    memset(powerSum, 0, sizeof(powerSum));
}

float Probe::power(int index) const {
    return powerSum[index] - ((index == 0) ? 0.0f : powerSum[index-1]);
}

vec3 Probe::sample(std::mt19937 &gen, const BRDF &brdf, float &prob) const {
    std::uniform_real_distribution<float> U(0.0f, 1.0f);
    auto sample_xy = [&](float &x, float &y, float &P) {
        float randomValue = U(gen) * powerSum[probeSize * probeSize - 1];
        int index = std::lower_bound(powerSum, powerSum + probeSize * probeSize, randomValue) - powerSum;
        int i = index / probeSize, j = index % probeSize;
        x = (i + U(gen)) / probeSize;
        y = (j + U(gen)) / probeSize;
        P = power(index) / Epower;
    };

    float x, y;
    for (int T = 0; T < maxTrys; T++) {
        float pdf_gird;
        sample_xy(x, y, pdf_gird);
        if (x*x + y*y > 1.0f)
            continue;
        vec3 dir = x * tangent + y * bitangent + sqrt(1.0f - x*x - y*y) * normal;
        float P_accept = brdf.P_accept(dir);
        if (P_accept > 0.0f && U(gen) < P_accept) {
            prob *= pdf_gird;
            return dir;
        }
    }
    return brdf.sample(gen);
}