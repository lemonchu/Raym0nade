#ifndef RAYM0NADE_PATH_MATERIAL_GLSL
#define RAYM0NADE_PATH_MATERIAL_GLSL

vec3 pathVertexNormal(uint vertexId) {
    const PackedVertex vertex = vertexBuffer.vertices[vertexId];
    return vec3(
        vertex.positionAndNormalX.w,
        vertex.normalYZAndUv.x,
        vertex.normalYZAndUv.y);
}

void pathMakeTangentSpace(
    vec3 normal, vec3 incoming, out vec3 tangent, out vec3 bitangent) {
    const vec3 unitNormal = safeNormalize(normal, vec3(0.0, 0.0, 1.0));
    const vec3 candidate = cross(unitNormal, incoming);
    if (!isFiniteVector(candidate) || dot(candidate, candidate) <= 1.0e-12) {
        const vec3 helper = abs(unitNormal.x) < 0.8
                                ? vec3(1.0, 0.0, 0.0)
                                : vec3(0.0, 1.0, 0.0);
        tangent = safeNormalize(cross(helper, unitNormal), vec3(0.0, 1.0, 0.0));
        bitangent = safeNormalize(cross(unitNormal, tangent), vec3(1.0, 0.0, 0.0));
        return;
    }
    bitangent = safeNormalize(candidate, vec3(1.0, 0.0, 0.0));
    tangent = safeNormalize(cross(bitangent, unitNormal), vec3(0.0, 1.0, 0.0));
}

PathRayDifferential pathInitialDifferential(vec3 unnormalizedDirection) {
    PathRayDifferential differential;
    differential.positionDx = vec3(0.0);
    differential.positionDy = vec3(0.0);
    differential.directionDx = vec3(0.0);
    differential.directionDy = vec3(0.0);

    const vec3 directionDx =
        pathParameters.cameraPositionAndPixelScale.w *
        pathParameters.cameraRightAndSampleWeight.xyz;
    const vec3 directionDy =
        pathParameters.cameraPositionAndPixelScale.w *
        pathParameters.cameraUpAndDirectProbability.xyz;
    const float squaredLength = dot(unnormalizedDirection, unnormalizedDirection);
    if (!(squaredLength > 1.0e-12) || !isFiniteScalar(squaredLength)) {
        return differential;
    }
    const float denominator = sqrt(squaredLength) * squaredLength;
    differential.directionDx =
        (squaredLength * directionDx -
         unnormalizedDirection * dot(unnormalizedDirection, directionDx)) /
        denominator;
    differential.directionDy =
        (squaredLength * directionDy -
         unnormalizedDirection * dot(unnormalizedDirection, directionDy)) /
        denominator;
    if (!isFiniteVector(differential.directionDx)) {
        differential.directionDx = vec3(0.0);
    }
    if (!isFiniteVector(differential.directionDy)) {
        differential.directionDy = vec3(0.0);
    }
    return differential;
}

void pathPositionDifferentials(
    vec3 rayDirection,
    float hitDistance,
    vec3 normal,
    PathRayDifferential base,
    out vec3 positionDx,
    out vec3 positionDy) {
    const float denominator = dot(rayDirection, normal);
    if (!isFiniteScalar(denominator) || abs(denominator) <= 1.0e-8) {
        positionDx = vec3(0.0);
        positionDy = vec3(0.0);
        return;
    }
    const float distanceDx =
        -dot(base.positionDx + hitDistance * base.directionDx, normal) /
        denominator;
    const float distanceDy =
        -dot(base.positionDy + hitDistance * base.directionDy, normal) /
        denominator;
    positionDx = base.positionDx + distanceDx * rayDirection +
                 hitDistance * base.directionDx;
    positionDy = base.positionDy + distanceDy * rayDirection +
                 hitDistance * base.directionDy;
    if (!isFiniteVector(positionDx)) {
        positionDx = vec3(0.0);
    }
    if (!isFiniteVector(positionDy)) {
        positionDy = vec3(0.0);
    }
}

PathRayDifferential pathReflectionDifferential(
    vec3 rayDirection,
    vec3 surfaceNormal,
    vec3 positionDx,
    vec3 positionDy,
    PathRayDifferential base) {
    PathRayDifferential result;
    result.positionDx = positionDx;
    result.positionDy = positionDy;
    result.directionDx =
        base.directionDx - 2.0 * dot(base.directionDx, surfaceNormal) * surfaceNormal;
    result.directionDy =
        base.directionDy - 2.0 * dot(base.directionDy, surfaceNormal) * surfaceNormal;
    if (!isFiniteVector(result.directionDx)) {
        result.directionDx = vec3(0.0);
    }
    if (!isFiniteVector(result.directionDy)) {
        result.directionDy = vec3(0.0);
    }
    return result;
}

float pathTextureFootprint(
    uint primitiveId, vec3 positionDx, vec3 positionDy) {
    const vec2 uvDx = textureDerivative(primitiveId, positionDx);
    const vec2 uvDy = textureDerivative(primitiveId, positionDy);
    if (!isFiniteVector(uvDx) || !isFiniteVector(uvDy)) {
        return 0.0;
    }
    const float footprint = 0.5 * (length(uvDx) + length(uvDy));
    return isFiniteScalar(footprint) && footprint > 0.0 ? footprint : 0.0;
}

vec3 pathEmissiveColor(
    PackedMaterial material, vec2 uv, float textureFootprint) {
    const vec3 factor = max(material.emissionAndIor.xyz, vec3(0.0));
    if ((material.flagsAndReserved.x & pathMaterialHasEmissiveTexture) == 0U) {
        return factor;
    }
    const uint textureId = material.textureIds.z;
    const vec3 encoded =
        samplePackedTextureAtFootprint(textureId, uv, textureFootprint).rgb;
    const vec3 linear = pow(max(encoded, vec3(0.0)), vec3(diffuseGamma));
    const vec3 emission = linear * factor;
    return isFiniteVector(emission) ? emission : vec3(0.0);
}

void pathSurfaceParameters(
    PackedMaterial material, vec2 uv, out float roughness, out float metallic) {
    roughness = material.transmissionAndRoughness.w;
    metallic = material.metallicSpecularAndReserved.x;
    if ((material.flagsAndReserved.x & pathMaterialHasSpecularTexture) != 0U) {
        const vec4 surfaceData = samplePackedMip(material.textureIds.y, uv, 0U);
        metallic = clamp(surfaceData.b, 0.0, 0.99);
        roughness = clamp(surfaceData.g, 1.0e-3, 1.0);
    }
    roughness = isFiniteScalar(roughness) ? clamp(roughness, 1.0e-3, 1.0) : 0.8;
    metallic = isFiniteScalar(metallic) ? clamp(metallic, 0.0, 1.0) : 0.0;
}

vec3 pathApplyNormalMap(
    uint primitiveId,
    vec3 shapeNormal,
    vec3 surfaceNormal,
    vec3 normalMap) {
    if (!isFiniteVector(normalMap) ||
        dot(normalMap, normalMap) <= pathRayEpsilon * pathRayEpsilon) {
        return surfaceNormal;
    }
    const uvec3 vertexIds = triangleVertexIds(primitiveId);
    const vec3 v0 = vertexPosition(vertexIds.x);
    const vec3 edge1 = vertexPosition(vertexIds.y) - v0;
    const vec3 edge2 = vertexPosition(vertexIds.z) - v0;
    const vec2 uv0 = vertexUv(vertexIds.x);
    const vec2 deltaUv1 = vertexUv(vertexIds.y) - uv0;
    const vec2 deltaUv2 = vertexUv(vertexIds.z) - uv0;
    const float determinant = deltaUv1.x * deltaUv2.y -
                              deltaUv2.x * deltaUv1.y;
    if (!isFiniteScalar(determinant) || abs(determinant) <= 1.0e-8) {
        return surfaceNormal;
    }
    const float inverseDeterminant = 1.0 / determinant;
    const vec3 tangentCandidate =
        inverseDeterminant * (deltaUv2.y * edge1 - deltaUv1.y * edge2);
    const vec3 tangent = normalizedOrZero(
        tangentCandidate - shapeNormal * dot(shapeNormal, tangentCandidate));
    if (dot(tangent, tangent) <= 0.0) {
        return surfaceNormal;
    }
    const vec3 bitangentCandidate =
        inverseDeterminant * (-deltaUv2.x * edge1 + deltaUv1.x * edge2);
    const float handedness =
        dot(cross(shapeNormal, tangent), bitangentCandidate) < 0.0 ? -1.0 : 1.0;
    const vec3 bitangent =
        handedness * safeNormalize(cross(shapeNormal, tangent), vec3(0.0));
    const vec3 mapped = tangent * normalMap.x + bitangent * normalMap.y +
                        surfaceNormal * normalMap.z;
    return safeNormalize(mapped, surfaceNormal);
}

bool pathPopulateSurface(
    PackedTraceHit hit,
    vec3 rayOrigin,
    vec3 rayDirection,
    PathRayDifferential differential,
    out vec3 positionDx,
    out vec3 positionDy,
    out PathSurface surface) {
    surface.position = rayOrigin + hit.distance * rayDirection;
    surface.baseColor = vec3(0.0);
    surface.emission = vec3(0.0);
    surface.specular = 0.04;
    surface.roughness = 0.8;
    surface.metallic = 0.0;
    surface.opacity = 1.0;
    surface.eta = 1.0;
    surface.materialId = triangleMaterialBuffer.materialIds[hit.primitiveId];
    surface.entering = true;

    const uvec3 vertexIds = triangleVertexIds(hit.primitiveId);
    const vec3 v0 = vertexPosition(vertexIds.x);
    const vec3 v1 = vertexPosition(vertexIds.y);
    const vec3 v2 = vertexPosition(vertexIds.z);
    surface.shapeNormal = safeNormalize(cross(v1 - v0, v2 - v0), -rayDirection);
    surface.entering = dot(surface.shapeNormal, -rayDirection) >= 0.0;
    if (!surface.entering) {
        surface.shapeNormal = -surface.shapeNormal;
    }
    if (dot(surface.shapeNormal, surface.shapeNormal) <= 0.0) {
        surface.surfaceNormal = surface.shapeNormal;
        positionDx = vec3(0.0);
        positionDy = vec3(0.0);
        return false;
    }

    vec3 n0 = safeNormalize(pathVertexNormal(vertexIds.x), surface.shapeNormal);
    vec3 n1 = safeNormalize(pathVertexNormal(vertexIds.y), surface.shapeNormal);
    vec3 n2 = safeNormalize(pathVertexNormal(vertexIds.z), surface.shapeNormal);
    if (dot(n0, surface.shapeNormal) < 0.0) n0 = -n0;
    if (dot(n1, surface.shapeNormal) < 0.0) n1 = -n1;
    if (dot(n2, surface.shapeNormal) < 0.0) n2 = -n2;
    const vec3 weights = barycentricWeights(hit.barycentrics);
    surface.surfaceNormal = safeNormalize(
        weights.x * n0 + weights.y * n1 + weights.z * n2,
        surface.shapeNormal);
    const vec3 originalSurfaceNormal = surface.surfaceNormal;

    pathPositionDifferentials(
        rayDirection,
        hit.distance,
        surface.shapeNormal,
        differential,
        positionDx,
        positionDy);
    const float footprint =
        pathTextureFootprint(hit.primitiveId, positionDx, positionDy);
    const vec2 uv = interpolateUv(hit.primitiveId, hit.barycentrics);
    const PackedMaterial material = materialBuffer.materials[surface.materialId];
    pathSurfaceParameters(material, uv, surface.roughness, surface.metallic);
    surface.opacity = material.diffuseAndOpacity.w;
    surface.eta = material.emissionAndIor.w;
    surface.specular = material.metallicSpecularAndReserved.y;
    if (surface.opacity > 1.0 - pathRayEpsilon) {
        surface.entering = true;
    }
    surface.baseColor = evaluateBaseColor(
        hit.primitiveId, hit.barycentrics, footprint);
    surface.emission = pathEmissiveColor(material, uv, footprint);

    if ((material.flagsAndReserved.x & pathMaterialHasNormalTexture) != 0U) {
        const uint textureId = material.textureIds.w;
        const vec3 normalMap =
            samplePackedTextureAtFootprint(textureId, uv, footprint).rgb *
            2.0 - 1.0;
        surface.surfaceNormal = pathApplyNormalMap(
            hit.primitiveId,
            surface.shapeNormal,
            surface.surfaceNormal,
            normalMap);
    }
    if (dot(surface.surfaceNormal, rayDirection) >= 0.0) {
        surface.surfaceNormal = originalSurfaceNormal;
    }
    if (dot(surface.surfaceNormal, rayDirection) >= 0.0) {
        surface.surfaceNormal = surface.shapeNormal;
    }

    return isFiniteVector(surface.position) &&
           isFiniteVector(surface.shapeNormal) &&
           isFiniteVector(surface.surfaceNormal) &&
           isFiniteVector(surface.baseColor) &&
           isFiniteVector(surface.emission) &&
           isFiniteScalar(surface.specular) &&
           isFiniteScalar(surface.roughness) &&
           isFiniteScalar(surface.metallic) &&
           isFiniteScalar(surface.opacity) &&
           isFiniteScalar(surface.eta);
}

#endif
