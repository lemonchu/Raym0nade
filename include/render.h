#ifndef RENDER_H
#define RENDER_H

#include <iostream>
#include "model.h"
#include "image.h"
#include "geometry.h"

#define MAX_RAY_DEPTH 16

struct RenderArgs {
    vec3 position, direction, up, right;
    float accuracy;     // 胶片距离为 1.0，每个像素的偏移量为 accuracy
    unsigned int oversampling, width, height;
    std::string savePath;
};

struct HitInfo {
    float t;
    glm::vec3 normal, emission;
    vec4 diffuseColor;
    HitInfo() : t(0), normal(vec3(0)), emission(vec3(0)), diffuseColor(vec4(0)) {}
};

class Renderer {
public:
    Model *modelPtr;
    Renderer();

    HitInfo Renderer::rayHit(Ray ray);
    PixelData sampleRay(Ray ray, int depth);

    void render(Model &model, const RenderArgs &args);
};

#endif // RENDER_H