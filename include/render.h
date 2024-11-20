#ifndef RENDER_H
#define RENDER_H

#include "model.h"
#include "image.h"
#include "geometry.h"

struct RenderArgs {
    vec3 position, direction, up, right;
    float accuracy;     // 胶片距离为 1.0，每个像素的偏移量为 accuracy
    unsigned int oversampling, width, height;
    std::string savePath;
};

class Renderer {
public:
    Model *modelPtr;
    Renderer();

    PixelData sampleRay(Ray ray, int depth);

    void render(Model &model, const RenderArgs &args);
};

#endif // RENDER_H