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

    PixelData sampleRay(Ray ray, int depth = 0) {
        HitRecord hit = model.kdt.rayHit(ray);
        PixelData ret;
        if (hit.t_max < INFINITY) {

            const auto& tri = *hit.tri;
            const auto& texture = *tri.texture;
            const auto& uv = tri.uv;

            const auto& intersection = ray.origin + hit.t_max * ray.direction;

            float area = length(cross(tri.v[1] - tri.v[0], tri.v[2] - tri.v[0]));
            float area1 = length(cross(intersection - tri.v[0], tri.v[2] - tri.v[0])) / area;
            float area2 = length(cross(tri.v[1] - tri.v[0], intersection - tri.v[0])) / area;
            float area3 = 1.0f - area1 - area2;

            float a = area1, b = area2, c = area3;

            float u = a * uv[0].x + b * uv[1].x + c * uv[2].x;
            float v = a * uv[0].y + b * uv[1].y + c * uv[2].y;

            int texX = static_cast<int>(u * texture.width);
            int texY = texture.height - static_cast<int>(v * texture.height);

            const std::vector<uint8_t>& imageData = texture.getImage(1);
            int pixelIndex = (texY * texture.width + texX) * 3; // 3 channels for RGB

            ret.color = glm::vec<3, float>(
                    imageData[pixelIndex] / 255.0f,
                    imageData[pixelIndex + 1] / 255.0f,
                    imageData[pixelIndex + 2] / 255.0f
            );

            ret.depth = hit.t_max;

            vec3 light = normalize(vec3(0, -1, -1));
            float lightIntensity = 0.6;
            float ambientIntensity = 0.4;
            vec3 normal = normalize(cross(tri.v[1] - tri.v[0], tri.v[2] - tri.v[0]));
            float diffuse = glm::dot(normal, light);

            ret.color = ret.color * (ambientIntensity + lightIntensity * std::max(0.0f, diffuse));
        } else {
            ret.color = glm::vec<3, float>(1, 1, 1);
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