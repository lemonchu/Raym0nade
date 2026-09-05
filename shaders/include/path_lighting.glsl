#ifndef RAYM0NADE_PATH_LIGHTING_GLSL
#define RAYM0NADE_PATH_LIGHTING_GLSL

bool pathHasEnvironment() {
    return pathParameters.seedLightCountAndEnvironmentWidthHeight.z > 0U &&
           pathParameters.seedLightCountAndEnvironmentWidthHeight.w > 0U;
}

vec3 pathEnvironmentRadiance(vec3 direction) {
    const uint width =
        pathParameters.seedLightCountAndEnvironmentWidthHeight.z;
    const uint height =
        pathParameters.seedLightCountAndEnvironmentWidthHeight.w;
    if (width == 0U || height == 0U || !isFiniteVector(direction)) {
        return vec3(0.0);
    }
    const vec3 unitDirection = normalizedOrZero(direction);
    if (dot(unitDirection, unitDirection) == 0.0) {
        return vec3(0.0);
    }
    float azimuth = atan(-unitDirection.x, unitDirection.z);
    if (azimuth < 0.0) {
        azimuth += 2.0 * pi;
    }
    const float polar = acos(clamp(unitDirection.y, -1.0, 1.0));
    const uint column = min(
        uint(azimuth * float(width) / (2.0 * pi)), width - 1U);
    const uint row = min(uint(polar * float(height) / pi), height - 1U);
    const vec3 radiance =
        environmentTexelBuffer.environmentTexels[row * width + column]
            .radianceAndConditionalCdf.xyz;
    return isFiniteVector(radiance) ? max(radiance, vec3(0.0)) : vec3(0.0);
}

uint pathSelectEnvironmentRow(float target, uint height) {
    uint first = 0U;
    uint count = height;
    while (count > 0U) {
        const uint step = count / 2U;
        const uint middle = first + step;
        const float cdf =
            environmentRowBuffer.environmentRows[middle]
                .probabilityCdfSolidAngleAndReserved.y;
        if (!(target < cdf)) {
            first = middle + 1U;
            count -= step + 1U;
        } else {
            count = step;
        }
    }
    return min(first, height - 1U);
}

uint pathSelectEnvironmentColumn(float target, uint row, uint width) {
    const uint rowOffset = row * width;
    uint first = 0U;
    uint count = width;
    while (count > 0U) {
        const uint step = count / 2U;
        const uint middle = first + step;
        const float cdf =
            environmentTexelBuffer.environmentTexels[rowOffset + middle]
                .radianceAndConditionalCdf.w;
        if (!(target < cdf)) {
            first = middle + 1U;
            count -= step + 1U;
        } else {
            count = step;
        }
    }
    return min(first, width - 1U);
}

bool pathSampleEnvironment(
    PathSurface surface,
    vec3 incoming,
    float categoryProbability,
    float selectionRandom,
    out vec3 throughput,
    out vec3 radiance) {
    throughput = vec3(0.0);
    radiance = vec3(0.0);
    const uint width =
        pathParameters.seedLightCountAndEnvironmentWidthHeight.z;
    const uint height =
        pathParameters.seedLightCountAndEnvironmentWidthHeight.w;
    const uint flags =
        pathParameters.environmentFlagsTotalSppSceneVersionAndFlags.x;
    if (width == 0U || height == 0U ||
        (flags & pathEnvironmentHasImportance) == 0U ||
        !(categoryProbability > 0.0)) {
        return false;
    }

    const uint row = pathSelectEnvironmentRow(selectionRandom, height);
    const vec4 rowData =
        environmentRowBuffer.environmentRows[row]
            .probabilityCdfSolidAngleAndReserved;
    const float previousRowCdf = row == 0U
                                     ? 0.0
                                     : environmentRowBuffer.environmentRows[row - 1U]
                                           .probabilityCdfSolidAngleAndReserved.y;
    const float rowProbability = rowData.y - previousRowCdf;
    if (!(rowProbability > 0.0) || !isFiniteScalar(rowProbability) ||
        !(rowData.z > 0.0) || !isFiniteScalar(rowData.z)) {
        return false;
    }
    const float conditionalTarget = clamp(
        (selectionRandom - previousRowCdf) / rowProbability,
        0.0,
        0.99999994);
    const uint column =
        pathSelectEnvironmentColumn(conditionalTarget, row, width);
    const uint texelIndex = row * width + column;
    const vec4 texel =
        environmentTexelBuffer.environmentTexels[texelIndex]
            .radianceAndConditionalCdf;
    const float previousColumnCdf = column == 0U
                                        ? 0.0
                                        : environmentTexelBuffer
                                              .environmentTexels[texelIndex - 1U]
                                              .radianceAndConditionalCdf.w;
    const float conditionalProbability = texel.w - previousColumnCdf;
    const float probability = rowProbability * conditionalProbability;
    if (!(probability > 0.0) || !isFiniteScalar(probability)) {
        return false;
    }

    const float polar = pi * (float(row) + 0.5) / float(height);
    const float azimuth = 2.0 * pi * (float(column) + 0.5) / float(width);
    const vec3 direction = vec3(
        -sin(polar) * sin(azimuth),
        cos(polar),
        sin(polar) * cos(azimuth));
    const PackedTraceHit shadowHit = tracePackedScene(
        surface.position,
        pathParameters.cameraDirectionAndRayMinimum.w,
        direction,
        rayMaximum,
        gl_RayFlagsTerminateOnFirstHitEXT);
    if (shadowHit.hasHit) {
        return false;
    }
    throughput = pathEvaluateBsdf(surface, incoming, direction);
    radiance = texel.xyz *
               (rowData.z / probability) /
               categoryProbability;
    return isFiniteVector(throughput) && isFiniteVector(radiance);
}

float pathAreaLightWeight(vec3 position, PackedAreaLight light) {
    const vec3 offset = light.centerAndPower.xyz - position;
    const float distanceSquared = dot(offset, offset);
    const float power = light.centerAndPower.w;
    if (!isFiniteScalar(distanceSquared) || !isFiniteScalar(power) ||
        !(power > 0.0)) {
        return 0.0;
    }
    return power / max(distanceSquared, pathMinimumLightDistanceSquared);
}

uint pathSelectAreaLight(
    vec3 position, float selectionRandom, out float objectProbability) {
    objectProbability = 0.0;
    const uint lightCount =
        pathParameters.seedLightCountAndEnvironmentWidthHeight.y;
    float totalWeight = 0.0;
    for (uint index = 0U; index < lightCount; ++index) {
        totalWeight += pathAreaLightWeight(
            position, areaLightBuffer.areaLights[index]);
    }
    if (!(totalWeight > 0.0) || !isFiniteScalar(totalWeight)) {
        return pathInvalidSceneId;
    }

    const float target = totalWeight * selectionRandom;
    float cumulative = 0.0;
    uint lastValid = pathInvalidSceneId;
    float selectedWeight = 0.0;
    for (uint index = 0U; index < lightCount; ++index) {
        const float weight = pathAreaLightWeight(
            position, areaLightBuffer.areaLights[index]);
        if (!(weight > 0.0)) {
            continue;
        }
        lastValid = index;
        cumulative += weight;
        if (target < cumulative) {
            selectedWeight = weight;
            objectProbability = selectedWeight / totalWeight;
            return index;
        }
    }
    if (lastValid != pathInvalidSceneId) {
        selectedWeight = pathAreaLightWeight(
            position, areaLightBuffer.areaLights[lastValid]);
        objectProbability = selectedWeight / totalWeight;
    }
    return lastValid;
}

uint pathSelectAreaTriangle(
    PackedAreaLight light, float selectionRandom) {
    const uint firstTriangle = light.triangleRangeAndReserved.x;
    const uint triangleCount = light.triangleRangeAndReserved.y;
    if (triangleCount == 0U) {
        return pathInvalidSceneId;
    }
    uint first = 0U;
    uint count = triangleCount;
    while (count > 0U) {
        const uint step = count / 2U;
        const uint middle = first + step;
        const float cdf =
            areaLightTriangleBuffer.areaLightTriangles[firstTriangle + middle]
                .areaProbabilityCdfAndReserved.z;
        if (!(selectionRandom < cdf)) {
            first = middle + 1U;
            count -= step + 1U;
        } else {
            count = step;
        }
    }
    return firstTriangle + min(first, triangleCount - 1U);
}

vec3 pathAreaTriangleEmission(
    PackedAreaLightTriangle triangle, vec3 barycentrics) {
    const vec2 uv =
        barycentrics.x * vertexUv(triangle.vertexIdsAndMaterialId.x) +
        barycentrics.y * vertexUv(triangle.vertexIdsAndMaterialId.y) +
        barycentrics.z * vertexUv(triangle.vertexIdsAndMaterialId.z);
    const PackedMaterial material =
        materialBuffer.materials[triangle.vertexIdsAndMaterialId.w];
    return pathEmissiveColor(material, uv, 0.0);
}

bool pathSampleAreaLight(
    PathSurface surface,
    vec3 incoming,
    float categoryProbability,
    float objectRandom,
    float triangleRandom,
    float firstRandom,
    float secondRandom,
    out vec3 throughput,
    out vec3 radiance) {
    throughput = vec3(0.0);
    radiance = vec3(0.0);
    if (!(categoryProbability > 0.0)) {
        return false;
    }
    float objectProbability = 0.0;
    const uint objectId = pathSelectAreaLight(
        surface.position, objectRandom, objectProbability);
    if (objectId == pathInvalidSceneId || !(objectProbability > 0.0)) {
        return false;
    }
    const PackedAreaLight light = areaLightBuffer.areaLights[objectId];
    const uint triangleId = pathSelectAreaTriangle(light, triangleRandom);
    if (triangleId == pathInvalidSceneId) {
        return false;
    }
    const PackedAreaLightTriangle triangle =
        areaLightTriangleBuffer.areaLightTriangles[triangleId];
    const float area = triangle.areaProbabilityCdfAndReserved.x;
    const uint firstTriangle = light.triangleRangeAndReserved.x;
    const float previousTriangleCdf = triangleId == firstTriangle
                                          ? 0.0
                                          : areaLightTriangleBuffer
                                                .areaLightTriangles[triangleId - 1U]
                                                .areaProbabilityCdfAndReserved.z;
    const float faceProbability =
        triangle.areaProbabilityCdfAndReserved.z - previousTriangleCdf;
    if (!(area > pathMinimumPdf) || !(faceProbability > 0.0) ||
        !isFiniteScalar(area) || !isFiniteScalar(faceProbability)) {
        return false;
    }
    const float root = pathSafeSqrt(firstRandom);
    const vec3 barycentrics = vec3(
        1.0 - root,
        root * (1.0 - secondRandom),
        root * secondRandom);
    const vec3 v0 = vertexPosition(triangle.vertexIdsAndMaterialId.x);
    const vec3 v1 = vertexPosition(triangle.vertexIdsAndMaterialId.y);
    const vec3 v2 = vertexPosition(triangle.vertexIdsAndMaterialId.z);
    const vec3 lightPosition =
        barycentrics.x * v0 + barycentrics.y * v1 + barycentrics.z * v2;
    const vec3 offset = lightPosition - surface.position;
    const float distanceSquared = dot(offset, offset);
    if (!(distanceSquared > pathMinimumLightDistanceSquared) ||
        !isFiniteScalar(distanceSquared)) {
        return false;
    }
    const float distance = pathSafeSqrt(distanceSquared);
    const vec3 direction = offset / distance;
    const vec3 lightNormal = safeNormalize(cross(v1 - v0, v2 - v0), vec3(0.0));
    const float lightCosine = max(0.0, dot(lightNormal, -direction));
    if (!(lightCosine > pathMinimumPdf)) {
        return false;
    }
    const float shadowMaximum =
        max(distance - pathRayEpsilon,
            pathParameters.cameraDirectionAndRayMinimum.w);
    const PackedTraceHit shadowHit = tracePackedScene(
        surface.position,
        pathParameters.cameraDirectionAndRayMinimum.w,
        direction,
        shadowMaximum,
        gl_RayFlagsTerminateOnFirstHitEXT);
    if (shadowHit.hasHit) {
        return false;
    }
    const float areaPdf = objectProbability * faceProbability / area;
    const float solidAnglePdf = areaPdf * distanceSquared / lightCosine;
    const float combinedPdf = categoryProbability * solidAnglePdf;
    if (!(combinedPdf > pathMinimumPdf) || !isFiniteScalar(combinedPdf)) {
        return false;
    }
    throughput = pathEvaluateBsdf(surface, incoming, direction);
    radiance = pathAreaTriangleEmission(triangle, barycentrics) / combinedPdf;
    return isFiniteVector(throughput) && isFiniteVector(radiance);
}

bool pathSampleDirectLight(
    PathSurface surface,
    vec3 incoming,
    uint pixelIndex,
    uint sampleIndex,
    uint bounceIndex,
    uint replicateIndex,
    uint lightSampleIndex,
    out vec3 throughput,
    out vec3 radiance) {
    throughput = vec3(0.0);
    radiance = vec3(0.0);
    const bool hasEnvironment = pathHasEnvironment();
    const bool hasAreaLights =
        pathParameters.seedLightCountAndEnvironmentWidthHeight.y > 0U;
    if (!hasEnvironment && !hasAreaLights) {
        return false;
    }
    const float environmentProbability =
        hasEnvironment && hasAreaLights ? 0.5 : (hasEnvironment ? 1.0 : 0.0);
    const float areaProbability = 1.0 - environmentProbability;
    const uint base = pathLightDimension(lightSampleIndex, 0U);
    const float categoryRandom = pathRandom(
        pixelIndex,
        sampleIndex,
        bounceIndex,
        replicateIndex,
        base + pathDimensionLightCategory);
    const bool chooseEnvironment = hasEnvironment &&
        (!hasAreaLights || categoryRandom < environmentProbability);
    bool valid = false;
    if (chooseEnvironment) {
        const float selectionRandom = pathRandom(
            pixelIndex,
            sampleIndex,
            bounceIndex,
            replicateIndex,
            base + pathDimensionLightObject);
        valid = pathSampleEnvironment(
            surface,
            incoming,
            environmentProbability,
            selectionRandom,
            throughput,
            radiance);
    } else {
        valid = pathSampleAreaLight(
            surface,
            incoming,
            areaProbability,
            pathRandom(
                pixelIndex,
                sampleIndex,
                bounceIndex,
                replicateIndex,
                base + pathDimensionLightObject),
            pathRandom(
                pixelIndex,
                sampleIndex,
                bounceIndex,
                replicateIndex,
                base + pathDimensionLightTriangle),
            pathRandom(
                pixelIndex,
                sampleIndex,
                bounceIndex,
                replicateIndex,
                base + pathDimensionLightFirst),
            pathRandom(
                pixelIndex,
                sampleIndex,
                bounceIndex,
                replicateIndex,
                base + pathDimensionLightSecond),
            throughput,
            radiance);
    }
    if (valid) {
        pathLimitThroughput(throughput);
    }
    return valid && isFiniteVector(throughput) && isFiniteVector(radiance);
}

#endif
