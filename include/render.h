#ifndef RENDER_H
#define RENDER_H

#include <glm/glm.hpp>
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
    Renderer() {}

    PixelData sampleRay(Ray ray, int depth = 0) {
        HitRecord hit = modelPtr->kdt.rayHit(ray);
        PixelData ret;
        if (hit.t_max < INFINITY) {

            const auto& tri = *hit.tri;
            const auto& texture = *tri.texture;

            /*if (texture.enabled[1]) {
                const auto& uv = tri.uv;

                const auto& intersection = ray.origin + hit.t_max * ray.direction;

                float
                        area = length(cross(tri.v[1] - tri.v[0], tri.v[2] - tri.v[0])),
                        area2 = length(cross(tri.v[1] - tri.v[0], intersection - tri.v[0])) / area,
                        area1 = length(cross(intersection - tri.v[0], tri.v[2] - tri.v[0])) / area,
                        area0 = 1.0f - area1 - area2;
                vec2 texUV = area0 * uv[0] + area1 * uv[1] + area2 * uv[2];

                int texX = lround(texUV[0] * texture.width + 0.5);
                int texY = texture.height - lround(texUV[1] * texture.height + 0.5);

                const std::vector<uint8_t>& imageData = texture.getImage(1);
                int pixelIndex = (texY * texture.width + texX) * 3; // 3 channels for RGB

                ret.color = glm::vec<3, float>(
                        imageData[pixelIndex] / 255.0f,
                        imageData[pixelIndex + 1] / 255.0f,
                        imageData[pixelIndex + 2] / 255.0f
                );

            } else {
                ret.color = glm::vec<3, float>(1, 1, 1);
            }*/

            ret.depth = hit.t_max;
            ret.color = glm::vec<3, float>(-log(ret.depth));


            /*vec3 light = normalize(vec3(0, -1, -1));
            float lightIntensity = 0.6;
            float ambientIntensity = 0.4;
            vec3 normal = normalize(cross(tri.v[1] - tri.v[0], tri.v[2] - tri.v[0]));
            float diffuse = glm::dot(normal, light);
            ret.color = ret.color * (ambientIntensity + lightIntensity * std::max(0.0f, diffuse));*/
        } else {
            ret.color = glm::vec<3, float>(1,0,0);
            ret.depth = 0;
        }

        return ret;
    }

    void render(Model &model, const RenderArgs &args) {
        modelPtr = &model;
        const unsigned int
            width = args.width,
            height = args.height,
            Oversampling = args.oversampling;
        const float accuracy = args.accuracy;
        vec3
            direction = args.direction,
            right = args.right,
            up = args.up,
            position = args.position;
        Image image(width, height);
        for (int x = 0; x < width; x++)
            for (int y = 0; y < height; y++) {
                PixelData &pixel = image.buffer[y * width + x];
                for (int x_os = 0; x_os < Oversampling; x_os++)
                    for (int y_os = 0; y_os < Oversampling; y_os++) {
                        float
                            rayX = x + 1.0 * x_os / Oversampling - width / 2,
                            rayY = y + 1.0 * y_os / Oversampling - height / 2;
                        glm::vec<3, float> aim =
                                direction + accuracy * (rayX * right + rayY * up);
                        aim = normalize(aim);
                        Ray ray = {position, aim};
                        pixel = pixel + sampleRay(ray);
                    }
                pixel.color /= Oversampling * Oversampling;
                pixel.depth /= Oversampling * Oversampling;
            }
        image.normalize();
        image.save(args.savePath.c_str());
    }
};

#endif // RENDER_H