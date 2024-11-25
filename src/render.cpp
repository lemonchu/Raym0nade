#include <random>
#include <corecrt_math_defines.h>
#include <iostream>
#include "render.h"

Renderer::Renderer() {}

glm::vec3 barycentric(const glm::vec3& A, const glm::vec3& B, const glm::vec3& C, const glm::vec3& P) {
    glm::vec3 v0 = B - A;
    glm::vec3 v1 = C - A;
    glm::vec3 v2 = P - A;
    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;
    return glm::vec3(u, v, w);
}

const float areaThereshold = 1e-2;

HitInfo getHitInfo(const Face& face, const glm::vec3& intersection) {
    HitInfo hitInfo;
    const glm::vec3 baryCoords = barycentric(face.v[0], face.v[1], face.v[2], intersection);
    const Material& material = *face.material;
    vec2 texUV =
            baryCoords[0] * face.data[0]->uv
            + baryCoords[1] * face.data[1]->uv
            + baryCoords[2] * face.data[2]->uv;
    if (material.isEnabled(TextureIdForDiffuseColor))
        hitInfo.diffuseColor = material.getImage(TextureIdForDiffuseColor).get(texUV[0], texUV[1]);
    vec3 cx = cross(face.v[1] - face.v[0], face.v[2] - face.v[0]);
    if (length(cx) < areaThereshold) {
        hitInfo.normal = normalize(
                baryCoords[0] * face.data[0]->normal
                + baryCoords[1] * face.data[1]->normal
                + baryCoords[2] * face.data[2]->normal
        );
    } else {
        hitInfo.normal = normalize(cx);
    }
    hitInfo.emission = face.material->emission;
    if (length(hitInfo.emission) > 0) {
        vec3 emissionColor = material.getImage(TextureIdForEmission).get(texUV[0], texUV[1]);
        hitInfo.emission *= emissionColor;
    }
    return hitInfo;
}

std::mt19937 gen;

HitInfo Renderer::rayHit(Ray ray) {
    HitRecord hit = modelPtr->kdt.rayHit(ray);
    if (hit.t_max < INFINITY) {
        const auto& face = *hit.face;
        HitInfo hitInfo = getHitInfo(face, ray.origin + ray.direction * hit.t_max);
        hitInfo.t = hit.t_max;
        return hitInfo;
    }
    return HitInfo();
}

/*
        // Apply lighting
        vec3 light = normalize(vec3(10, -1, -1));
        float specularIntensity = 0.7;
        float diffuseIntensity = 0.7;

        // Calculate diffuse term
        float diffuseFactor = std::max(dot(normal, light), 0.0f);
        vec3 diffuseTerm = diffuseIntensity * diffuseFactor * diffuseColor;

        // Calculate specular term
        vec3 viewDir = normalize(-ray.direction);
        vec3 reflectDir = reflect(-light, normal);
        float spec = pow(std::max(dot(viewDir, reflectDir), 0.0f), shininess);
        vec3 specularTerm = specularIntensity * spec * specularColor;

        // Final color calculation
        ret.color = diffuseColor;//diffuseTerm + specularTerm;
        // std::cout << "Color: " << ret.color.x << " " << ret.color.y << " " << ret.color.z << std::endl;*/

//ret.depth = hit.t_max;
//ret.color = - vec3(log(ret.depth)); */

const int maxRandomRounds = 16;
const float eps_lightRadius = 1e-3;
const unsigned long long minTouchSamples = 2048;

vec3 sampleDirectLight(const vec3 &intersection, const HitInfo &info, Model &model) {
    auto &lightObjects = model.lightObjects;

    /*vec3 light = vec3(0.0f);
    for (const auto& lightObject : lightObjects) {
        float distance = glm::length(lightObject.center - intersection);
        float dot = std::max(glm::dot(normalize(lightObject.center - intersection), info.normal), 0.00f);
        light += dot * lightObject.color * lightObject.power / (distance * distance + eps_lightRadius);
    }
    return light * 10.0f;*/

    float totalWeight = 0.0f;
    std::vector<float> weights;
    for (const auto& lightObject : lightObjects) {
        float distance = glm::length(lightObject.center - intersection);
        float dot = std::max(glm::dot(normalize(lightObject.center - intersection), info.normal), 0.05f);
        float weight = dot * lightObject.power / (distance * distance * distance + eps_lightRadius);
        if (lightObject.total >= minTouchSamples)
            weight *= 1.0f * lightObject.touch / lightObject.total;
        weights.push_back(weight);
        totalWeight += weight;
    }

    if (totalWeight == 0.0f)
        return vec3(-1.0f);

    std::discrete_distribution<int> lightDist(weights.begin(), weights.end());
    int lightIndex = lightDist(gen);
    auto& lightObject = lightObjects[lightIndex];
    int faceIndex = lightObject.faceDist(gen);
    auto& lightFace = lightObject.lightFaces[faceIndex];
    vec3 lightDir;
    float cosTheta, cosPhi;
    for (int T = 0; T < maxRandomRounds; T++) {
        lightDir = normalize(lightFace.position - intersection);
        cosTheta = glm::dot(info.normal, lightDir);
        cosPhi = glm::dot(lightFace.normal, -lightDir);
        if (cosTheta > eps_zero && cosPhi > eps_zero)
            break;
        faceIndex = lightObject.faceDist(gen);
        lightFace = lightObject.lightFaces[faceIndex];
    }
    if (cosTheta < eps_zero || cosPhi < eps_zero)
        return vec3(-1.0f);
    float distance = glm::length(lightFace.position - intersection);

    Ray shadowRay = {intersection + eps_zero * lightDir, lightDir};
    HitRecord shadowHit = model.kdt.rayHit(shadowRay);

    lightObject.total++;

    if (shadowHit.t_max < distance - eps_lightRadius)
        return vec3(0.0f);

    lightObject.touch++;

    float distanceSquared = distance * distance + eps_lightRadius;
    float lightObjectProbability = weights[lightIndex] / totalWeight;
    float lightFaceProbability = lightFace.power / lightObject.power;
    float pdf = lightObjectProbability * lightFaceProbability;
    float contribution = lightFace.power * cosPhi * cosTheta / (distanceSquared * 2.0f * M_PI * pdf);
    if (rand()%20000==0)
        std::cout << "Contribution: " << contribution << std::endl;
    return lightObject.color * contribution * 100.0f;
    return contribution == 0.0 ? vec3(0.0f) : vec3(1.0);
}

const float transparentDepth = 1e-6;

vec3 Renderer::sampleRay(Ray ray, int depth = 0) {
    HitInfo info = rayHit(ray);
    vec3 intersection = ray.origin + ray.direction * info.t;

    if (info.t == INFINITY)
        return vec3(0.5);

    if (depth < MAX_RAY_DEPTH && info.diffuseColor[3] == 0.0f) {
        Ray ray2 = {ray.origin + ray.direction * (info.t + transparentDepth), ray.direction};
        return sampleRay(ray2, depth + 1);
    }

    if (length(info.emission) > 0)
        return vec3(1, 0, 0);
    else {
        vec3 color = vec3(info.diffuseColor[0], info.diffuseColor[1], info.diffuseColor[2]);
        vec3 directLight = sampleDirectLight(intersection, info, *modelPtr);
        if (directLight == vec3(-1.0f))
            return vec3(-1.0f);
        color *= sampleDirectLight(intersection, info, *modelPtr);
        return color;
    }
}

void Renderer::render(Model &model, const RenderArgs &args) {
    modelPtr = &model;
    const unsigned int
        width = args.width,
        height = args.height,
        spp = args.spp,
        oversampling = args.oversampling;
    const float
        accuracy = args.accuracy;
    vec3
        direction = args.direction,
        right = args.right,
        up = args.up,
        position = args.position;
    Image image(width, height);
    for (unsigned int T = 0; T < spp; T++) {
        for (unsigned int x = 0; x < width; x++)
            for (unsigned int y = 0; y < height; y++) {
                for (unsigned int x_os = 0; x_os < oversampling; x_os++)
                    for (unsigned int y_os = 0; y_os < oversampling; y_os++) {
                        float
                                rayX = x + 1.0f * x_os / oversampling - width / 2.0f,
                                rayY = y + 1.0f * y_os / oversampling - height / 2.0f;
                        vec3 aim =
                                direction + accuracy * (rayX * right + rayY * up);
                        aim = normalize(aim);
                        Ray ray = {position, aim};
                        vec3 color = sampleRay(ray);
                        if (color != vec3(-1.0f)) {
                            image.buffer[y * width + x].color += color;
                            image.buffer[y * width + x].cnt ++;
                        }
                    }
            }
    }
    image.save(args.savePath.c_str());

    for (auto &lightObject : model.lightObjects) {
        std::cout << "Light object " << &lightObject - &model.lightObjects[0] << " touched " << lightObject.touch << " times, total " << lightObject.total << " times." << std::endl;
        std::cout << "Touch probability: " << 1.0f * lightObject.touch / lightObject.total << std::endl;
    }
}