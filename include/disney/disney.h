/*
 * Disney BSDF Shading Model
 *
 * This file is part of the Selas rendering engine.
 * Original repository: https://github.com/schuttejoe/Selas
 *
 * Modifications have been made to this file to enhance the BSDF implementation.
 *
 * Copyright (c) 2018 Joe Schutte
 *
 * Modifications Copyright (c) 2024 Meng Chu and Yuchen Yang
 *
 * The original code is licensed under the MIT License. See the LICENSE file or
 * visit the original repository for details.
 */

#ifndef RAYM0NADE_DISNEY_H
#define RAYM0NADE_DISNEY_H

#include <glm/glm.hpp>
#include <cmath>

#include "sampling.h"
#include "disney/Fresnel.h"
#include "disney/utils.h"
#include "component.h"

class Generator;
struct BsdfSample;

// -- BSDF evaluation for next event estimation
glm::vec3 EvaluateDisney(const BRDF &surface, glm::vec3 outDir, bool thin, float &forwardPdf, float &reversePdf);

// -- Shaders
bool SampleDisney(Generator &gen, const BRDF &surface, glm::vec3 v, bool thin, BsdfSample &sample);

#endif //RAYM0NADE_DISNEY_H
