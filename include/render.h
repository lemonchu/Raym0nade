#ifndef RENDER_H
#define RENDER_H

#include "model.h"
#include "geometry.h"

struct RenderArgs {
    vec3 position, direction, up, right;
    float accuracy, exposure; // 胶片距离为 1.0，每个像素的偏移量为 accuracy
    unsigned int width, height, spp, threads;
    std::string savePath;
};

struct RenderData {
    std::mt19937 gen;
    int T_RayAndTexture;
    explicit RenderData(int seed);
};

void render_multiThread(Model &model, const RenderArgs &args);

#endif // RENDER_H