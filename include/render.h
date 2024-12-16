#ifndef RENDER_H
#define RENDER_H

#include "model.h"
#include "geometry.h"

struct RenderArgs {
    vec3 position, direction, up, right;
    float accuracy, exposure, P_Direct; // 胶片距离为 1.0，每个像素的偏移量为 accuracy
    int width, height, spp, threads;
    std::string savePath;
};

struct RenderData {
    Generator gen;
    int C_lightSamples;
    std::atomic<int> *renderedPixels;
    explicit RenderData(int seed);
};

void render_multiThread(Model &model, const RenderArgs &args);

#endif // RENDER_H