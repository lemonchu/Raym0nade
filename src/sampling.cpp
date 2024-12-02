#include "sampling.h"
#include "sobel.h"

const float M_PI = 3.14159265359f;

BRDF::BRDF(const vec3 &inDir, const vec3 &shapeNormal, const vec3 &surfaceNormal) :
        inDir(inDir), shapeNormal(shapeNormal), surfaceNormal(surfaceNormal), roughness(1.0f), max(1.0f) {
    getTangentSpaceWithInDir(surfaceNormal, inDir, tangent, bitangent);
    EP_accept = 0.5;
}

float BRDF::pdf(const vec3 &outDir) const {
    return dot(outDir,surfaceNormal) <= 0.0f ? 0.0f : std::max(dot(outDir, surfaceNormal), 0.0f);
}

float BRDF::P_accept(const vec3 &outDir) const {
    return pdf(outDir) / max;
}

vec3 BRDF::sample(Generator genU, Generator genV, float &P_success) const {
    auto sample0 = [&]() -> vec3 {
        float u = genU(), v = genV() * 2.0f * M_PI;
        float d = sqrt(u);
        float z = sqrt(1 - d*d);
        float x = d * cos(v), y = d * sin(v);
        return x * tangent + y * bitangent + z * surfaceNormal;
    };
    vec3 direction;
    for (int T = 1; T <= MaxTrys; T++) {
        direction = sample0();
        if (dot(direction, shapeNormal) > 0.0f) {
            P_success = 1.0f / T;
            return direction;
        }
    }
    P_success = 0.0f;
    return vec3(0.0f);
}

int sample(const std::vector<float> &weights, float randomValue) {
    float cumulativeWeight = 0.0f;
    for (int i = 0; i < weights.size(); ++i) {
        cumulativeWeight += weights[i];
        if (randomValue <= cumulativeWeight)
            return i;
    }
}

void sampleLightObject(const vec3 &pos, const BRDF &brdf, const Model &model, Generator &gen, float &prob, int &lightIndex) {
    float totalWeight = 0.0f;
    std::vector<float> weights;
    for (int id = 0; id < model.lightObjects.size(); id++)  {
        const LightObject &lightObject = model.lightObjects[id];
        vec3 lightDir = normalize(lightObject.center - pos);
        float distance = glm::length(lightObject.center - pos);
        float brdf_pdf = std::max(brdf.pdf(lightDir) + 0.1f, 0.0f);
        float weight = brdf_pdf * lightObject.power / (distance * distance + eps_lightRadius);
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

void sampleLightFace(const vec3 &pos, const LightObject &lightObject, Generator &genF, Generator &genP, float &P_success, int &faceIndex) {
    vec3 lightDir;
    for(int T = 1; T <= MaxTrys; T++) {
        faceIndex = lightObject.faceDist(genF);
        auto &lightFace = lightObject.lightFaces[faceIndex];
        lightDir = normalize(lightFace.position - pos);
        float cosPhi = dot(lightFace.normal, -lightDir);
        if (cosPhi <= 0.0f || genP() > cosPhi) {
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

vec3 sampleDirectLight(const vec3 &pos, const BRDF &brdf, const Model &model, Generator genU, Generator genV, Generator genP) {

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
    sampleLightObject(pos, brdf,  model, genU, P_lightObject, lightIndex);
    if (lightIndex == -1)
        return vec3(0.0f);
    auto& lightObject = model.lightObjects[lightIndex];

    float P_success;
    int faceIndex;
    sampleLightFace(pos, lightObject, genV, genP, P_success, faceIndex);
    if (faceIndex == -1)
        return vec3(0.0f);
    auto &lightFace = lightObject.lightFaces[faceIndex];

    vec3 lightDir = normalize(lightFace.position - pos);
    float distance = glm::length(lightFace.position - pos);
    if (model.rayHit_test({pos, lightDir}, distance - eps_lightRadius))
        return vec3(0.0f);

    float distanceSquared = distance * distance + eps_lightRadius;
    float contribution = lightObject.power * brdf.pdf(lightDir) / (distanceSquared * P_lightObject) * P_success;
    vec3 light = lightObject.color * contribution;
    return light;
}
