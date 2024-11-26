#include <corecrt_math_defines.h>
#include <thread>
#include <random>
#include <iostream>
#include "render.h"

const float areaThereshold = 5e-3; // 小面平滑插值，大面直接取法向量

bool TransparentTest(const Ray &ray, const HitRecord &hit) {
    const Face& face = *hit.face;
    const Material& material = *face.material;
    if (!material.hasTransparentPart)
        return false;
    const vec3 intersection = ray.origin + ray.direction * hit.t_max;
    const vec3 baryCoords = barycentric(face.v[0], face.v[1], face.v[2], intersection);
    vec2 texUV =
            baryCoords[0] * face.data[0]->uv
            + baryCoords[1] * face.data[1]->uv
            + baryCoords[2] * face.data[2]->uv;
    vec4 diffuseColor = material.getDiffuseColor(texUV[0], texUV[1]);
    return diffuseColor[3] < 0.5f;
}

void getHitInfo(const Face& face, const glm::vec3& intersection, HitInfo &hitInfo) {
    const glm::vec3 baryCoords = barycentric(face.v[0], face.v[1], face.v[2], intersection);
    const Material& material = *face.material;
    vec2 texUV =
            baryCoords[0] * face.data[0]->uv
            + baryCoords[1] * face.data[1]->uv
            + baryCoords[2] * face.data[2]->uv;
    hitInfo.diffuseColor = material.getDiffuseColor(texUV[0], texUV[1]);
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
    if (face.lightObject != nullptr) {
        hitInfo.emission = face.lightObject->color * face.lightObject->powerDensity;
    } else
        hitInfo.emission = vec3(0.0f);
}

const int maxRayDepth_hit = 8;

bool rayHit_test(Ray ray, const Model &model, float aimDepth) {
    HitRecord hit(eps_zero, aimDepth + eps_zero);
    for (int T = 0; T < maxRayDepth_hit; T++) {
        model.kdt.rayHit(ray, hit);
        if (hit.t_max > aimDepth)
            return false;
        if (!TransparentTest(ray, hit))
            return true;
        hit = HitRecord(hit.t_max + eps_zero, aimDepth + eps_zero);
    }
    return true;
}

HitInfo rayHit(Ray ray, const Model &model) {
    HitRecord hit;
    HitInfo hitInfo;
    for (int T = 0; T < maxRayDepth_hit; T++) {
        model.kdt.rayHit(ray, hit);
        if (hit.t_max == INFINITY) {
            hitInfo.t = INFINITY;
            break;
        }
        vec3 intersection = ray.origin + ray.direction * hit.t_max;
        getHitInfo(*hit.face, intersection, hitInfo);
        if (hitInfo.diffuseColor[3] > 0.0f){
            hitInfo.t = hit.t_max;
            return hitInfo;
        }
        hit = HitRecord(hit.t_max + eps_zero, INFINITY);
    }
    return hitInfo;
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
const float Ecos = 0.5f;

vec4 sampleDirectLight(const vec3 &pos, const HitInfo &info, const Model &model, std::mt19937 &gen) {

    /*vec3 light = vec3(0.0f);
    auto &lightObjects = model.lightObjects;
    for (const auto& lightObject : lightObjects) {
        float distance = glm::length(lightObject.center - pos);
        float dot = std::max(glm::dot(normalize(lightObject.center - pos), info.normal), 0.00f);
        light += dot * lightObject.color * lightObject.power / (distance * distance + eps_lightRadius);
    }
    return vec4(light, 1.0f);*/

    float totalWeight = 0.0f;
    std::vector<float> weights;
    for (const auto& lightObject : model.lightObjects) {
        float distance = glm::length(lightObject.center - pos);
        float dot = std::max(glm::dot(normalize(lightObject.center - pos), info.normal), 0.05f);
        float weight = dot * lightObject.power / (distance * distance * distance + eps_lightRadius);
        weights.push_back(weight);
        totalWeight += weight;
    }
    if (totalWeight == 0.0f)
        return vec4(0.0f, 0.0f, 0.0f, 1.0f);
    std::discrete_distribution<int> lightDist(weights.begin(), weights.end());
    int lightIndex = lightDist(gen);
    auto& lightObject = model.lightObjects[lightIndex];

    int total = 0, faceIndex;
    vec3 lightDir;
    float cosTheta, cosPhi;
    for (int T = 0; T < maxRandomRounds; T++) {
        faceIndex = lightObject.faceDist(gen);
        auto &lightFace = lightObject.lightFaces[faceIndex];
        lightDir = normalize(lightFace.position - pos);
        cosTheta = glm::dot(info.normal, lightDir);
        cosPhi = glm::dot(lightFace.normal, -lightDir);
        if (cosTheta > eps_zero && cosPhi > eps_zero)
            break;
        else
            total++;
    }
    if (cosTheta < eps_zero || cosPhi < eps_zero)
        return vec4(0.0f, 0.0f, 0.0f, total);

    auto &lightFace = lightObject.lightFaces[faceIndex];
    float distance = glm::length(lightFace.position - pos);

    if (rayHit_test({pos, lightDir}, model, distance - eps_lightRadius))
        return vec4(0.0f, 0.0f, 0.0f, total+1);

    float distanceSquared = distance * distance + eps_lightRadius;
    float lightObjectProbability = weights[lightIndex] / totalWeight;
    float lightFaceProbability = lightFace.power / lightObject.power;
    float pdf = lightObjectProbability * lightFaceProbability;
    float contribution = lightFace.power * cosPhi * cosTheta / (distanceSquared * M_PI * pdf);
    return vec4(lightObject.color * contribution, total+1);
}

vec3 randomRayDirection(const vec3 &normal, std::mt19937 &gen) {
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float u = dist(gen), v = dist(gen) * 2.0f * M_PI;
    float d = sqrt(u);
    float z = sqrt(1 - d*d);
    float x = d * cos(v), y = d * sin(v);
    vec3 direction = vec3(x, y, z);
    tangentTransform(normal, direction);
    return direction;
}

const float StopProb = 0.36f;
const int maxRayDepth = 16;

vec4 sampleRay(Ray ray, std::mt19937 &gen, const Model &model) {
    vec4 ret = vec4(1.0f);
    float prob = 1.0f;
    for (int T = 0;; T++) {
        HitInfo info = rayHit(ray, model);

        if (info.t == INFINITY)
            return vec4(0.0f, 0.0f, 0.0f, 1.0f);

        if (length(info.emission) > 0) {
            if (T > 0)
                return vec4(0.0f, 0.0f, 0.0f, 1.0f);
            else
                return vec4(info.emission, 1.0f);
        }

        vec3 intersection = ray.origin + ray.direction * info.t;
        ret *= info.diffuseColor;

        if (T == maxRayDepth || std::uniform_real_distribution<float>(0.0f, 1.0f)(gen) < StopProb) {
            vec4 directLight = sampleDirectLight(intersection, info, model, gen);
            ret = ret * directLight;
            if (T < maxRayDepth)
                prob *= StopProb;
            ret.x /= prob;
            ret.y /= prob;
            ret.z /= prob;
            return ret;
        }

        vec3 newDirection = randomRayDirection(info.normal, gen);
        ray = {intersection, newDirection};
        prob *= (1.0f - StopProb);
    }
}

void render(const Model &model, const RenderArgs &args, std::mt19937 &gen, Image &image) {
    const unsigned int
            width = args.width,
            height = args.height,
            spp = args.spp,
            oversampling = args.oversampling;
    const float
            accuracy = args.accuracy,
            exposure = args.exposure;
    vec3
            direction = args.direction,
            right = args.right,
            up = args.up,
            position = args.position;
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
                        vec4 color = sampleRay(ray, gen, model);
                        color[3] /= exposure;
                        image.buffer[y * width + x].color += color;
                    }
            }
        std::cout << "Frame " << T << " completed." << std::endl;
    }
}

void render_multithread(Model &model, const RenderArgs &args) {
    const unsigned int
            startTime = clock(),
            threads = args.threads;
    std::cout << "Rendering started with " << threads << " threads." << std::endl;

    std::vector<Image> images;
    std::vector<std::mt19937> gens;
    images.reserve(threads);
    gens.reserve(threads);
    for (unsigned int i = 0; i < threads; ++i) {
        images.emplace_back(args.width, args.height);
        gens.emplace_back(i);
    }
    std::vector<std::thread> threadPool;

    auto renderTask = [&](int threadIndex) {
        render(model, args, gens[threadIndex], images[threadIndex]);
    };

    for (unsigned int i = 0; i < threads; ++i) {
        threadPool.emplace_back(renderTask, i);
    }

    for (auto &thread : threadPool) {
        thread.join();
    }

    Image finalImage(args.width, args.height);
    for (const auto &image : images)
        for (unsigned int id = 0; id < args.height * args.width; ++id)
            finalImage.buffer[id].color += image.buffer[id].color;

    finalImage.save(args.savePath.c_str());

    std::cout << "Rendering completed in " << clock() - startTime << " ms." << std::endl;
}