#include "render.h"

Renderer::Renderer() {}

const int TextureIdForDiffuseColor = 1;
const int TextureIdForSpecularColor = 2;
const int TextureIdForAmbientColor = 3;
const int TextureIdForEmission = 4;

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


HitInfo getHitInfo(const Face& face, const glm::vec3& intersection) {
    HitInfo hitInfo;
    const glm::vec3 baryCoords = barycentric(face.v[0], face.v[1], face.v[2], intersection);
    const Material& material = *face.material;
    vec2 texUV =
            baryCoords[0] * face.data[0]->uv
            + baryCoords[1] * face.data[1]->uv
            + baryCoords[2] * face.data[2]->uv;
    if (material.isEnabled(TextureIdForDiffuseColor))
        hitInfo.diffuseColor = material.getImage(TextureIdForDiffuseColor).get(texUV[0], texUV[1]);
    hitInfo.normal = normalize(
            baryCoords[0] * face.data[0]->normal
            + baryCoords[1] * face.data[1]->normal
            + baryCoords[2] * face.data[2]->normal
    );
    hitInfo.emission = face.material->emission;
    if (length(hitInfo.emission) > 0) {
        vec3 emissionColor = material.getImage(TextureIdForEmission).get(texUV[0], texUV[1]);
        hitInfo.emission *= emissionColor;
    }
    return hitInfo;
}

HitInfo Renderer::rayHit(Ray ray) {
    HitRecord hit = modelPtr->kdt.rayHit(ray);
    if (hit.t_max < INFINITY) {
        const auto& face = *hit.face;
        HitInfo hitInfo = getHitInfo(face, ray.origin + ray.direction * hit.t_max);
        hitInfo.t = hit.t_max;
        return hitInfo;
    }
    return HitInfo();
}

/*
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
        ret.color = diffuseColor;//diffuseTerm + specularTerm;
        // std::cout << "Color: " << ret.color.x << " " << ret.color.y << " " << ret.color.z << std::endl;*/

//ret.depth = hit.t_max;
//ret.color = - vec3(log(ret.depth)); */

PixelData Renderer::sampleRay(Ray ray, int depth = 0) {
    HitInfo info = rayHit(ray);
    PixelData ret;
    ret.depth = info.t;
    if (length(info.emission) > 0)
        ret.color = vec3(1, 0, 0);
    else
        ret.color = vec3(info.diffuseColor[0], info.diffuseColor[1], info.diffuseColor[2]);
    if (depth < MAX_RAY_DEPTH && info.diffuseColor[3] < 1.0 - eps_zero) {
        Ray ray2 = {ray.origin + ray.direction * (info.t + eps_zero), ray.direction};
        PixelData next = sampleRay(ray2, depth + 1);
        ret.color = info.diffuseColor[3]*ret.color + (1-info.diffuseColor[3])*next.color;
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
    //image.normalize();
    image.save(args.savePath.c_str());
}