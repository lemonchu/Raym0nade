#include "render.h"

Renderer::Renderer() {}

const int TextureIdForDiffuseColor = 1;

glm::vec3 barycentric(const glm::vec3& A, const glm::vec3& B, const glm::vec3& C, const glm::vec3& P) {
    glm::vec3 v0 = B - A;
    glm::vec3 v1 = C - A;
    glm::vec3 v2 = P - A;
    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;
    return glm::vec3(u, v, w);
}

PixelData Renderer::sampleRay(Ray ray, int depth = 0) {
    HitRecord hit = modelPtr->kdt.rayHit(ray);
    PixelData ret;
    if (hit.t_max < INFINITY) {
        const auto& face = *hit.face;
        const auto intersection = ray.origin + hit.t_max * ray.direction;
        glm::vec3 baryCoords = barycentric(face.v[0], face.v[1], face.v[2], intersection);
        float u = baryCoords.x, v = baryCoords.y, w = baryCoords.z;

        const auto normal = normalize(u * face.data[0]->normal + v * face.data[1]->normal + w * face.data[2]->normal);
        ret.depth = hit.t_max;
        ret.color = vec3(0.0f, 0.0f, 0.0f); // Initialize default color to black

        const Material& texture = *face.texture;

        // Retrieve material properties
        vec3 diffuseColor = texture.diffuseColor;
        vec3 specularColor = texture.specularColor;
        vec3 ambientColor = texture.ambientColor;
        float shininess = texture.shininess;

        // Update diffuse color with texture if available
        if (texture.enabled[TextureIdForDiffuseColor]) {
            vec2 texUV = u * face.data[0]->uv + v * face.data[1]->uv + w * face.data[2]->uv;
            int texX = lround(texUV[0] * texture.width + 0.5);
            int texY = texture.height - lround(texUV[1] * texture.height + 0.5);
            texX = std::max(0, std::min(texX, texture.width - 1));
            texY = std::max(0, std::min(texY, texture.height - 1));

            const std::vector<uint8_t>& imageData = texture.getImage(1);
            int pixelIndex = (texY * texture.width + texX) * 3; // 3 channels for RGB
            diffuseColor = vec3(
                    imageData[pixelIndex] / 255.0f,
                    imageData[pixelIndex + 1] / 255.0f,
                    imageData[pixelIndex + 2] / 255.0f
            );
        }

        // Apply lighting
        vec3 light = normalize(vec3(10, -1, -1));
        float specularIntensity = 0.7;
        float diffuseIntensity = 0.7;

        // Calculate diffuse term
        float diffuseFactor = std::max(dot(normal, light), 0.0f);
        vec3 diffuseTerm = diffuseIntensity * diffuseFactor * diffuseColor;

        // Calculate specular term
        vec3 viewDir = normalize(-ray.direction);
        vec3 reflectDir = reflect(-light, normal);
        float spec = pow(std::max(dot(viewDir, reflectDir), 0.0f), shininess);
        vec3 specularTerm = specularIntensity * spec * specularColor;

        // Final color calculation
        ret.color = diffuseTerm + specularTerm;
        // std::cout << "Color: " << ret.color.x << " " << ret.color.y << " " << ret.color.z << std::endl;
    } else {
        ret.color = vec3(1, 0, 0); // Background color (red)
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