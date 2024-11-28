#ifndef RENDER_H
#define RENDER_H

#include "model.h"
#include "geometry.h"

struct RenderArgs {
    vec3 position, direction, up, right;
    float accuracy, exposure;     // 胶片距离为 1.0，每个像素的偏移量为 accuracy
    unsigned int width, height, oversampling, spp, threads, probes;
    std::string savePath;
};

void render_multiThread(Model &model, const RenderArgs &args);

#endif // RENDER_H