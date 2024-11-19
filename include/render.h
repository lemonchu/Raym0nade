#ifndef RENDER_H
#define RENDER_H

#include <glm/glm.hpp>
#include "model.h"
#include "image.h"

struct Camera {
    glm::vec<3, float> position, direction, up, right;
    float accuracy;     // 胶片距离为 1.0，每个像素的偏移量为 accuracy
    int Oversampling;   // 超采样倍数
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

    void render() {
        for (int x = 0; x < image.width; x++)
            for (int y = 0; y < image.height; y++) {
                PixelData &pixel = image.buffer[y * image.width + x];
                for (int x_os = 0; x_os < camera.Oversampling; x_os++)
                    for (int y_os = 0; y_os < camera.Oversampling; y_os++) {
                        float
                            rayX = x + 1.0 * x_os / camera.Oversampling - image.width / 2,
                            rayY = y + 1.0 * y_os / camera.Oversampling - image.height / 2;
                        glm::vec<3, float> aim =
                                camera.direction + camera.accuracy * (rayX * camera.right + rayY * camera.up);
                        aim = normalize(aim);
                        Ray ray = {camera.position, aim};
                        pixel = pixel + sampleRay(ray);
                    }
                pixel.color /= camera.Oversampling * camera.Oversampling;
                pixel.depth /= camera.Oversampling * camera.Oversampling;
            }
        image.normalize();
    }

    void clearImage() {
        image.clear();
    }

    void saveImage(const char* file_name) {
        image.save(file_name);
    }
};

#endif // RENDER_H