#include <thread>
#include <random>
#include <iostream>
#include "render.h"
#include "octree.h"
#include "image.h"

const float P_RR = 0.7f;
const int maxRayDepth = 16;

vec4 sampleRay(Ray ray, const Model &model, const Octree &octree, std::mt19937 &gen) {
    vec3 radiance = vec4(1.0f);
    float prob = 1.0f;
    int failCount = 0;
    for (int T = 0;; T++) {
        HitInfo info = model.rayHit(ray);

        if (info.t == INFINITY)
            return vec4(0.0f, 0.0f, 0.0f, failCount+1);

        if (length(info.emission) > 0)
            return (T==0) ?
                vec4(info.emission, failCount+1) :
                vec4(0.0f, 0.0f, 0.0f, failCount+1);

        vec3 intersection = ray.origin + ray.direction * info.t;
        radiance *= (vec3)info.diffuseColor;
        BRDF brdf(ray.direction, info.shapeNormal, info.surfaceNormal);

        if (T == maxRayDepth || std::uniform_real_distribution<float>(0.0f, 1.0f)(gen) > P_RR) {
            vec3 light;
            sampleDirectLight(intersection, brdf, model, gen, failCount, light);
            radiance *= light;
            if (T < maxRayDepth)
                prob *= (1.0f - P_RR);
            return vec4(radiance / prob, failCount+1);
        }

        vec3 newDirection = brdf.sample(gen);
        ray = {intersection, newDirection};
        prob *= P_RR;
    }
}

void shootImportron(Ray ray, const Model &model, std::mt19937 &gen, std::vector<Importron> &importrons) {
    vec4 ret = vec4(1.0f);
    float weight = 1.0f;
    for (int T = 0;; T++) {
        HitInfo info = model.rayHit(ray);

        if (info.t == INFINITY)
            return ;

        vec3 intersection = ray.origin + ray.direction * info.t;
        float diffuseRate = (info.diffuseColor[0] + info.diffuseColor[1] + info.diffuseColor[2]) / 3.0f;
        if (T == maxRayDepth || std::uniform_real_distribution<float>(0.0f, 1.0f)(gen) > P_RR) {
            if (T < maxRayDepth)
                weight /= (1.0f - P_RR);
            importrons.push_back({intersection, info.surfaceNormal, weight});
            return ;
        }
        BRDF brdf(ray.direction, info.shapeNormal, info.surfaceNormal);
        vec3 newDirection = brdf.sample(gen);
        ray = {intersection, newDirection};
        weight *= diffuseRate / P_RR;
    }
}

const int minDepth = 6, maxDepth = 10;
const float coordLimit = 32.0f;

void BuildOctree(const Model &model, const RenderArgs &args, Octree &octree) {
    const unsigned int
            width = args.width,
            height = args.height,
            probes = args.probes;
    const float
            accuracy = args.accuracy;
    vec3
            direction = args.direction,
            right = args.right,
            up = args.up,
            position = args.position;
    std::vector<Importron> importrons;
    importrons.reserve(width * height);
    std::mt19937 gen(0);
    for (unsigned int x = 0; x < width; x++)
        for (unsigned int y = 0; y < height; y++) {
            float
                    rayX = x - width / 2.0f,
                    rayY = y - height / 2.0f;
            vec3 aim =
                    normalize(direction + accuracy * (rayX * right + rayY * up));
            Ray ray = {position, aim};
            shootImportron(ray, model, gen, importrons);
        }
    float totalWeight = 0.0f;
    for (const Importron &importron : importrons)
        totalWeight += importron.weight;

    std::cout << "Importrons: " << importrons.size() << std::endl;
    std::cout << "Total weight: " << totalWeight << std::endl;

    octree.coordLimit = coordLimit;
    octree.minWeight = totalWeight / probes * (9.0f/16);
    octree.minDepth = minDepth;
    octree.maxDepth = maxDepth;
    octree.build(importrons);

    std::cout << "Octree has built." << std::endl;
}

void render(const Model &model, const Octree &octree, const RenderArgs &args, std::mt19937 &gen, Image &image, int xL, int xR) {
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
        for (unsigned int x = xL; x < xR; x++)
            for (unsigned int y = 0; y < height; y++) {
                for (unsigned int x_os = 0; x_os < oversampling; x_os++)
                    for (unsigned int y_os = 0; y_os < oversampling; y_os++) {
                        float
                                rayX = x + 1.0f * x_os / oversampling - width / 2.0f,
                                rayY = y + 1.0f * y_os / oversampling - height / 2.0f;
                        vec3 aim =
                                normalize(direction + accuracy * (rayX * right + rayY * up));
                        Ray ray = {position, aim};
                        vec4 color = sampleRay(ray, model, octree, gen);
                        color[3] /= exposure;
                        image.buffer[y * width + x].color += color;
                    }
            }
        std::cout << "Frame " << T << " completed." << std::endl;
    }
}

void render_multiThread(Model &model, const RenderArgs &args) {
    const unsigned int
            startTime = clock(),
            width = args.width,
            threads = args.threads;

    Octree octree;
    /*BuildOctree(model, args, octree);*/

    std::cout << "Rendering started with " << threads << " threads." << std::endl;

    Image image(width, args.height);
    std::vector<std::mt19937> gens;
    gens.reserve(threads);
    for (unsigned int i = 0; i < threads; ++i)
        gens.emplace_back(0);
    std::vector<std::thread> threadPool;
    auto renderTask = [&](int threadIndex) {
        render(model, octree, args, gens[threadIndex], image,
               width/threads * threadIndex, width/threads * (threadIndex+1));
    };
    for (unsigned int i = 0; i < threads; ++i)
        threadPool.emplace_back(renderTask, i);
    for (auto &thread : threadPool)
        thread.join();

    image.save(args.savePath.c_str());

    std::cout << "Rendering completed in " << clock() - startTime << " ms." << std::endl;
}