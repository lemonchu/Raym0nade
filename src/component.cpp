#include "component.h"

Generator::Generator(unsigned int seed) :
    mt(seed), U(0.0f, 1.0f) {}

float Generator::operator()() {
    return U(mt);
}

void RandomDistribution::Init(const std::vector<float>& distribution) {
    prefixSums.resize(distribution.size());
    prefixSums[0] = distribution[0];
    for (size_t i = 1; i < distribution.size(); ++i) {
        prefixSums[i] = prefixSums[i - 1] + distribution[i];
    }
}

int RandomDistribution::operator()(Generator &gen) const {
    float randomValue = prefixSums.back() * gen();
    return std::lower_bound(prefixSums.begin(), prefixSums.end(), randomValue) - prefixSums.begin();
}

LightObject::LightObject() : center(glm::vec3(0)), color(glm::vec3(0)), power(0) {}

VertexData::VertexData() = default;

VertexData::VertexData(const vec2 &uv, const vec3 &normal) : uv(uv), normal(normal) {}

vec3 Face::center() const {
    return (v[0] + v[1] + v[2]) / 3.0f;
}

Box Face::aabb() const {
    return Box(
            vec3(
                    std::min(v[0][0],std::min(v[1][0],v[2][0])),
                    std::min(v[0][1],std::min(v[1][1],v[2][1])),
                    std::min(v[0][2],std::min(v[1][2],v[2][2]))
            ),
            vec3(
                    std::max(v[0][0],std::max(v[1][0],v[2][0])),
                    std::max(v[0][1],std::max(v[1][1],v[2][1])),
                    std::max(v[0][2],std::max(v[1][2],v[2][2]))
            )
    );
}

HitRecord::HitRecord() : t_min(eps_zero), t_max(INFINITY), face(nullptr) {}

HitRecord::HitRecord(float t_min, float t_max) : t_min(t_min), t_max(t_max), face(nullptr) {}