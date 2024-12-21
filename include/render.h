#ifndef RENDER_H
#define RENDER_H

#include "model.h"
#include "geometry.h"
#include <map>

struct RenderArgs {
    vec3 position, direction, up, right;
    float accuracy, exposure, P_Direct; // 胶片距离为 1.0，每个像素的偏移量为 accuracy
    int width, height, spp, threads;
    std::string savePath;
};

struct MediumData {
    float ior;
    vec3 absorb;
    MediumData(float ior, const vec3 &absorb);
};

class Medium {
private:
    std::multimap<int, MediumData> mediums;
public:
    void Init();
    void insert(int id, float ior, const vec3 &absorb);
    void erase(int id);
    [[nodiscard]] vec3 absorb() const;
    [[nodiscard]] float ior() const;
    [[nodiscard]] int size() const;
};

class RenderData {
public:
    Generator gen;
    Medium mediums;
    int C_lightSamples;
    std::atomic<int> *renderedPixels;
    explicit RenderData(int seed);
};

void render_multiThread(Model &model, const RenderArgs &args);

#endif // RENDER_H