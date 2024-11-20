#include "render.h"

Renderer::Renderer() {}

const int TextureIdForColor = 1;

vec3 getColor(const vec3& intersection, const Face& face) {
    float
        area = length(cross(face.v[1] - face.v[0], face.v[2] - face.v[0])),
        area2 = length(cross(face.v[1] - face.v[0], intersection - face.v[0])) / area,
        area1 = length(cross(intersection - face.v[0], face.v[2] - face.v[0])) / area,
        area0 = 1.0f - area1 - area2;
    const Texture& texture = *face.texture;
    if (texture.enabled[TextureIdForColor]) {
        vec2 texUV = area0 * face.data[0]->uv + area1 * face.data[1]->uv + area2 * face.data[2]->uv;

        int texX = lround(texUV[0] * texture.width + 0.5);
        int texY = texture.height - lround(texUV[1] * texture.height + 0.5);

        texX = std::max(0, std::min(texX, texture.width - 1));
        texY = std::max(0, std::min(texY, texture.height - 1));

        const std::vector<uint8_t>& imageData = texture.getImage(1);
        int pixelIndex = (texY * texture.width + texX) * 3; // 3 channels for RGB

        return vec3(
                imageData[pixelIndex] / 255.0f,
                imageData[pixelIndex + 1] / 255.0f,
                imageData[pixelIndex + 2] / 255.0f
        );
    } else {
        return area0 * face.data[0]->color + area1 * face.data[1]->color + area2 * face.data[2]->color;
    }
}


PixelData Renderer::sampleRay(Ray ray, int depth = 0) {
    HitRecord hit = modelPtr->kdt.rayHit(ray);
    PixelData ret;
    if (hit.t_max < INFINITY) {

        const auto& face = *hit.face;
        const auto& intersection = ray.origin + hit.t_max * ray.direction;

        ret.depth = hit.t_max;
        ret.color = getColor(intersection, face);

        /*vec3 light = normalize(vec3(0, -1, -1));
        float lightIntensity = 0.6;
        float ambientIntensity = 0.4;
        vec3 normal = normalize(cross(tri.v[1] - tri.v[0], tri.v[2] - tri.v[0]));
        float diffuse = glm::dot(normal, light);
        ret.color = ret.color * (ambientIntensity + lightIntensity * std::max(0.0f, diffuse));*/
    } else {
        ret.color = vec3(1,0,0);
        ret.depth = 0;
    }

    return ret;
}

void Renderer::render(Model &model, const RenderArgs &args) {
    modelPtr = &model;
    const unsigned int
        width = args.width,
        height = args.height,
        oversampling = args.oversampling;
    const float
        accuracy = args.accuracy,
        oversampleCnt = oversampling * oversampling;
    vec3
        direction = args.direction,
        right = args.right,
        up = args.up,
        position = args.position;
    Image image(width, height);
    for (unsigned int x = 0; x < width; x++)
        for (unsigned int y = 0; y < height; y++) {
            PixelData &pixel = image.buffer[y * width + x];
            for (unsigned int x_os = 0; x_os < oversampling; x_os++)
                for (unsigned int y_os = 0; y_os < oversampling; y_os++) {
                    float
                        rayX = x + 1.0f * x_os / oversampling - width / 2.0f,
                        rayY = y + 1.0f * y_os / oversampling - height / 2.0f;
                    vec3 aim =
                            direction + accuracy * (rayX * right + rayY * up);
                    aim = normalize(aim);
                    Ray ray = {position, aim};
                    pixel = pixel + sampleRay(ray);
                }
            pixel.color /= oversampleCnt;
            pixel.depth /= oversampleCnt;
        }
    // image.normalize();
    image.save(args.savePath.c_str());
}