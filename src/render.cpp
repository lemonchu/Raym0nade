#include <thread>
#include <random>
#include <iostream>
#include "render.h"

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

HitInfo getHitInfo(const Face& face, const glm::vec3& intersection) {
    HitInfo hitInfo;
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
    hitInfo.emission = face.material->emission;
    if (length(hitInfo.emission) > 0) {
        vec3 emissionColor = material.getImage(TextureIdForEmission).get(texUV[0], texUV[1]);
        hitInfo.emission *= emissionColor;
    }
    return hitInfo;
}

const int maxRayDepth = 8;

const float eps_lightRadius = 1e-3;

bool rayHit_test(Ray ray, const Model &model, float aimDepth) {
    float t = 0.0;
    for (int T = 0; T < maxRayDepth; T++) {
        HitRecord hit = model.kdt.rayHit(ray, aimDepth + eps_zero);
        t += hit.t_max;
        if (t > aimDepth)
            return false;
        if (TransparentTest(ray, hit))
            ray = {ray.origin + ray.direction * hit.t_max, ray.direction};
        else
            return true;
    }
    return true;
}

HitInfo rayHit(Ray ray, const Model &model) {
    float t = 0.0f;
    for (int T = 0; T < maxRayDepth; T++) {
        HitRecord hit = model.kdt.rayHit(ray);
        if (hit.t_max == INFINITY)
            break;
        HitInfo hitInfo = getHitInfo(*hit.face, ray.origin + ray.direction * hit.t_max);
        t += hit.t_max;
        if (hitInfo.diffuseColor[3] == 0.0f)
            ray = {ray.origin + ray.direction * (hit.t_max + eps_lightRadius), ray.direction};
        else {
            hitInfo.t = t;
            return hitInfo;
        }
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

vec4 sampleDirectLight(const vec3 &intersection, const HitInfo &info, const Model &model, std::mt19937 &gen) {
    auto &lightObjects = model.lightObjects;

    /*vec3 light = vec3(0.0f);
    for (const auto& lightObject : lightObjects) {
        float distance = glm::length(lightObject.center - intersection);
        float dot = std::max(glm::dot(normalize(lightObject.center - intersection), info.normal), 0.00f);
        light += dot * lightObject.color * lightObject.power / (distance * distance + eps_lightRadius);
    }
    return vec4(light, 1.0f);*/

    float totalWeight = 0.0f;
    std::vector<float> weights;
    for (const auto& lightObject : lightObjects) {
        float distance = glm::length(lightObject.center - intersection);
        float dot = std::max(glm::dot(normalize(lightObject.center - intersection), info.normal), 0.05f);
        float weight = dot * lightObject.power / (distance * distance * distance + eps_lightRadius);
        weights.push_back(weight);
        totalWeight += weight;
    }
    if (totalWeight == 0.0f)
        return vec4(0.0f, 0.0f, 0.0f, 1.0f);

    int total = 0;
    std::discrete_distribution<int> lightDist(weights.begin(), weights.end());
    int lightIndex = lightDist(gen);
    auto& lightObject = lightObjects[lightIndex];

    int faceIndex;
    vec3 lightDir;
    float cosTheta, cosPhi;
    for (int T = 0; T < maxRandomRounds; T++) {
        faceIndex = lightObject.faceDist(gen);
        auto &lightFace = lightObject.lightFaces[faceIndex];
        lightDir = normalize(lightFace.position - intersection);
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
    float distance = glm::length(lightFace.position - intersection);

    Ray shadowRay = {intersection + eps_zero * lightDir, lightDir};
    bool shadow = rayHit_test(shadowRay, model, distance - eps_lightRadius);

    if (shadow)
        return vec4(0.0f, 0.0f, 0.0f, total+1);

    float distanceSquared = distance * distance + eps_lightRadius;
    float lightObjectProbability = weights[lightIndex] / totalWeight;
    float lightFaceProbability = lightFace.power / lightObject.power;
    float pdf = lightObjectProbability * lightFaceProbability;
    float contribution = lightFace.power * cosPhi * cosTheta / (distanceSquared * pdf);
    return vec4(lightObject.color * contribution, total+1);
}

vec4 sampleRay(Ray ray, std::mt19937 &gen, const Model &model, int depth = 0) {
    HitInfo info = rayHit(ray, model);
    vec3 intersection = ray.origin + ray.direction * info.t;

    if (info.t == INFINITY)
        return vec4(0.0f, 0.0f, 1.0f, 1.0f);

    if (length(info.emission) > 0)
        return vec4(2.0f, 2.0f, 2.0f, 1.0f);
    else {
        vec4 directLight = sampleDirectLight(intersection, info, model, gen);
        return info.diffuseColor * directLight;
    }
}

void _render(const Model &model, const RenderArgs &args, std::mt19937 &gen, Image &image) {
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
        _render(model, args, gens[threadIndex], images[threadIndex]);
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