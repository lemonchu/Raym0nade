#ifndef MODEL_H
#define MODEL_H

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "material.h"
#include "component.h"
#include "kdt.h"
#include <Python.h>

struct HitInfo {
    vec3 shapeNormal, surfaceNormal, emission, baseColor, position;
    float specular, roughness, metallic, opacity, eta;
    int id;
    bool entering;
    HitInfo();
};

void getHitNormals(const Face& face, const vec3 &inDir, const vec3 &baryCoords,
                   vec3 &shapeNormal, vec3 &surfaceNormal_raw, bool &entering);
void getHitMaterial(const Face& face, const vec3 &baryCoords, const vec3 &hit_dPdx, const vec3 &hit_dPdy,
                    HitInfo &hitInfo);

class Model {
private:
    void processMaterial(const std::string &model_folder, const aiScene *scene);
    void checkLightObject(Face *meshFaces, aiMesh *mesh, const Material &material);
    void processMesh(aiMesh *mesh);
public:
    std::vector<Material> materials;
    std::vector<Face> faces;
    std::vector<VertexData> vertexDatas;
    std::vector<LightObject> lightObjects;
    KDT kdt;
    std::string model_path;

    Model();
    Model(const std::string &model_folder, const std::string &model_name);
    [[nodiscard]] HitRecord rayHit(Ray ray) const;
    [[nodiscard]] bool rayHit_test(Ray ray, float aimDepth) const;
};

#endif // MODEL_H