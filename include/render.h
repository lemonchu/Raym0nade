#ifndef RENDER_H
#define RENDER_H

#include "model.h"
#include "image.h"
#include "geometry.h"

#define MAX_RAY_DEPTH 16

struct RenderArgs {
    vec3 position, direction, up, right;
    float accuracy, exposure;     // 胶片距离为 1.0，每个像素的偏移量为 accuracy
    unsigned int width, height, oversampling, spp, threads;
    std::string savePath;
};

struct HitInfo {
    float t;
    glm::vec3 normal, emission;
    vec4 diffuseColor;
    HitInfo() : t(0), normal(vec3(0)), emission(vec3(0)), diffuseColor(vec4(0)) {}
};

void render_multithread(Model &model, const RenderArgs &args);

#endif // RENDER_H