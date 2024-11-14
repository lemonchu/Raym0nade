#ifndef RENDER_H
#define RENDER_H

#include <glm/glm.hpp>
#include "model.h"
#include "image.h"

struct Camera {
    glm::vec<3, float> position, direction, up, right;
    float accuracy;  // 胶片距离为 1.0，每个像素的偏移量为 accuracy
};

class Renderer {
private:
    Model model;
    Image image;
public:
    Camera camera;
    Renderer(int width, int height, const char* model_file) : image(width, height) {
        model.load(model_file);
    }
    HitRecord rayHit(Ray ray) {
        HitRecord closest_hit = HitRecord(vec3(0), eps, INFINITY);

        for (const Triangle &tri : model.triangles)
            RayTriangleIntersection(ray, tri, closest_hit);

        return closest_hit;
    }

    PixelData sampleRay(Ray ray, int depth = 0) {
        HitRecord hit = rayHit(ray);
        PixelData ret;
        if (hit.t_max < INFINITY) {
            ret.color = glm::vec<3, float>(1, 1, 1);
            ret.depth = hit.t_max;
        } else {
            ret.color = glm::vec<3, float>(0, 0, 0);
            ret.depth = 0;
        }
        return ret;
    }

    void render() {
        for (int x = 0; x < image.width; x++)
            for (int y = 0; y < image.height; y++) {
                glm::vec<3, float> aim =
                        + camera.direction
                        + camera.accuracy * (x-image.width/2) * camera.right
                        + camera.accuracy * (y-image.width/2) * camera.up;
                aim = normalize(aim);
                Ray ray = {camera.position, aim};
                image.buffer[y * image.width + x] = sampleRay(ray);
            }
    }

    void saveImage(const char* file_name) {
        image.save(file_name);
    }
};

#endif // RENDER_H