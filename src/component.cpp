#include <Python.h>
#include <iostream>
#include "component.h"

Generator::Generator(unsigned int seed) :
    mt(seed), U(1e-6f, 1.0f-1e-6f) {}

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
    return static_cast<int>
        (std::lower_bound(prefixSums.begin(), prefixSums.end(), randomValue) - prefixSums.begin());
}

float RandomDistribution::pdf(int index) const {
    float now = prefixSums[index];
    if (index > 0)
        now -= prefixSums[index - 1];
    return now / prefixSums.back();
}

LightObject::LightObject() : center(vec3(0)), color(vec3(0)), power(0) {}

VertexData::VertexData(const vec2 &uv, const vec3 &normal) : uv(uv), normal(normal) {}

vec3 Face::center() const {
    return (v[0] + v[1] + v[2]) / 3.0f;
}

Box Face::aabb() const {
    return Box(
            glm::min(v[0], glm::min(v[1], v[2])),
            glm::max(v[0], glm::max(v[1], v[2]))
    );
}

HitRecord::HitRecord() : t_min(eps_zero), t_max(INFINITY), face(nullptr) {}

HitRecord::HitRecord(float t_min, float t_max) : t_min(t_min), t_max(t_max), face(nullptr) {}

SkyBox::SkyBox() : width(0), height(0) {}

void SkyBox::Init() {
    std::vector<float> weights;
    weights.reserve(width * height);
    for (int v = 0; v < height; v++)
        for (int u = 0; u < width; u++) {
            float phi = PI * (float(v)+0.5f) / float(height);
            float area = sin(phi) * 2.0f * PI / float(width*height);
            int id = v * width + u;
            data[id] *= area;
            float Clum = dot(data[id], RGB_Weight);
            weights.emplace_back(Clum);
        }
    dist.Init(weights);
}

void SkyBox::load(const std::string &filename) {
    std::cout << "Loading HDR image from file: " << filename << std::endl;

    PyObject *pValue = call_image_to_array(filename, "hdr_to_array");

    if (pValue != nullptr) {
        PyObject *pWidth = PyTuple_GetItem(pValue, 0);
        PyObject *pHeight = PyTuple_GetItem(pValue, 1);
        PyObject *pChannels = PyTuple_GetItem(pValue, 2);
        PyObject *pData = PyTuple_GetItem(pValue, 3);

        if (pWidth == Py_None) {
            std::cerr << "Invalid width in Python object." << std::endl;
        }
        if (pHeight == Py_None) {
            std::cerr << "Invalid height in Python object." << std::endl;
        }
        if (pData == Py_None) {
            std::cerr << "Invalid data in Python object." << std::endl;
        }
        if (pWidth != Py_None && pHeight != Py_None && pData != Py_None) {
            width = PyLong_AsLong(pWidth);
            height = PyLong_AsLong(pHeight);
            std::cout << "Image dimensions: " << width << "x" << height << std::endl;

            Py_buffer view;
            if (PyObject_GetBuffer(pData, &view, PyBUF_SIMPLE) == 0) {
                float *buffer = static_cast<float*>(view.buf);
                int length = view.len / sizeof(vec3);
                data.reserve(length);
                for (int i = 0; i < length; i++)
                    data.emplace_back(buffer[i*3], buffer[i*3 + 1], buffer[i*3 + 2]);
                Init();
                PyBuffer_Release(&view);
            } else {
                std::cerr << "Failed to get buffer view from Python object." << std::endl;
            }
        } else {
            std::cerr << "Invalid data in Python object." << std::endl;
        }
        Py_DECREF(pValue);
    } else {
        PyErr_Print();
        std::cerr << "Failed to call Python function hdr_to_array." << std::endl;
    }
}

bool SkyBox::empty() const {
    return data.empty();
}

vec3 SkyBox::get(vec3 dir) const {
    if (data.empty())
        return vec3(0.0f);

    float theta = atan2(-dir.x, dir.z);
    float phi = acos(dir.y);
    if (theta < 0.0f) theta += 2.0f * PI;

    int u = static_cast<int>(theta / (2.0f * PI) * float(width));
    int v = static_cast<int>(phi / PI * float(height));

    if (u < 0) u = 0;
    if (u >= width) u = width - 1;
    if (v < 0) v = 0;
    if (v >= height) v = height - 1;

    phi = PI * (float(v)+0.5f) / float(height);
    float area = sin(phi) * 2.0f * PI / float(width*height);

    return data[v * width + u] / area;
}