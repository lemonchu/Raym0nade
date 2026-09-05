#ifndef RAYM0NADE_PATH_INTEGRATOR_GLSL
#define RAYM0NADE_PATH_INTEGRATOR_GLSL

const uint pathLightSampleCounts[17] = uint[](
    0U, 1U, 2U, 2U, 3U, 3U, 3U, 4U, 4U,
    4U, 4U, 5U, 5U, 5U, 5U, 5U, 6U);

void pathResetOutput(uint outputIndex) {
    pathOutputBuffer.pixels[outputIndex].shapeNormalAndRoughness = vec4(0.0);
    pathOutputBuffer.pixels[outputIndex].surfaceNormalAndSpecular = vec4(0.0);
    pathOutputBuffer.pixels[outputIndex].positionAndMetallic = vec4(0.0);
    pathOutputBuffer.pixels[outputIndex].baseColorAndOpacity = vec4(0.0);
    pathOutputBuffer.pixels[outputIndex].emissionAndEta = vec4(0.0);
    pathOutputBuffer.pixels[outputIndex]
        .materialEnteringHitAndDirectSampleCount =
        uvec4(pathInvalidSceneId, 1U, 0U, 0U);
    pathOutputBuffer.pixels[outputIndex].directDiffuse = vec4(0.0);
    pathOutputBuffer.pixels[outputIndex].directSpecular = vec4(0.0);
    pathOutputBuffer.pixels[outputIndex].indirectDiffuse = vec4(0.0);
    pathOutputBuffer.pixels[outputIndex].indirectSpecular = vec4(0.0);
}

void pathStoreSurface(uint outputIndex, PathSurface surface) {
    pathOutputBuffer.pixels[outputIndex].shapeNormalAndRoughness =
        vec4(surface.shapeNormal, surface.roughness);
    pathOutputBuffer.pixels[outputIndex].surfaceNormalAndSpecular =
        vec4(surface.surfaceNormal, surface.specular);
    pathOutputBuffer.pixels[outputIndex].positionAndMetallic =
        vec4(surface.position, surface.metallic);
    pathOutputBuffer.pixels[outputIndex].baseColorAndOpacity =
        vec4(surface.baseColor, surface.opacity);
    pathOutputBuffer.pixels[outputIndex].emissionAndEta =
        vec4(surface.emission, surface.eta);
    pathOutputBuffer.pixels[outputIndex]
        .materialEnteringHitAndDirectSampleCount =
        uvec4(
            surface.materialId,
            surface.entering ? 1U : 0U,
            1U,
            0U);
}

PathSurface pathLoadSurface(uint outputIndex) {
    const PathOutputPixel pixel = pathOutputBuffer.pixels[outputIndex];
    PathSurface surface;
    surface.shapeNormal = pixel.shapeNormalAndRoughness.xyz;
    surface.roughness = pixel.shapeNormalAndRoughness.w;
    surface.surfaceNormal = pixel.surfaceNormalAndSpecular.xyz;
    surface.specular = pixel.surfaceNormalAndSpecular.w;
    surface.position = pixel.positionAndMetallic.xyz;
    surface.metallic = pixel.positionAndMetallic.w;
    surface.baseColor = pixel.baseColorAndOpacity.xyz;
    surface.opacity = pixel.baseColorAndOpacity.w;
    surface.emission = pixel.emissionAndEta.xyz;
    surface.eta = pixel.emissionAndEta.w;
    surface.materialId = pixel.materialEnteringHitAndDirectSampleCount.x;
    surface.entering = pixel.materialEnteringHitAndDirectSampleCount.y != 0U;
    return surface;
}

void pathAccumulateRadiance(
    vec3 incomingRadiance, float weight, inout vec4 destination) {
    if (!isFiniteVector(incomingRadiance) || !isFiniteScalar(weight) ||
        !(weight > 0.0)) {
        return;
    }
    const vec3 weightedRadiance = incomingRadiance * weight;
    const vec3 accumulatedRadiance = destination.xyz + weightedRadiance;
    const float squaredMagnitude = dot(incomingRadiance, incomingRadiance);
    const float accumulatedSecondMoment =
        destination.w + squaredMagnitude * weight;
    if (!isFiniteVector(weightedRadiance) ||
        !isFiniteVector(accumulatedRadiance) ||
        !isFiniteScalar(accumulatedSecondMoment)) {
        return;
    }
    destination = vec4(accumulatedRadiance, accumulatedSecondMoment);
}

void pathAccumulateInwardRadiance(
    vec3 baseColor,
    vec3 sampleRadiance,
    vec3 sampleThroughput,
    float sampleWeight,
    inout vec4 diffuseRadiance,
    inout vec4 specularRadiance) {
    if (!isFiniteVector(baseColor) || !isFiniteVector(sampleRadiance) ||
        !isFiniteVector(sampleThroughput) || !isFiniteScalar(sampleWeight) ||
        !(sampleWeight > 0.0) || length(sampleRadiance) < pathRayEpsilon) {
        return;
    }
    const float baseLength = length(baseColor);
    if (!isFiniteScalar(baseLength) || baseLength < pathRayEpsilon) {
        pathAccumulateRadiance(
            sampleRadiance * sampleThroughput,
            sampleWeight,
            specularRadiance);
        return;
    }

    const float inverseSqrtThree = 0.5773502691896258;
    const vec3 white = vec3(inverseSqrtThree);
    const vec3 normalizedBase = safeNormalize(baseColor, vec3(0.0));
    const float alignment = clamp(dot(normalizedBase, white), -1.0, 1.0);
    const float minimumComponent = min(
        abs(baseColor.x), min(abs(baseColor.y), abs(baseColor.z)));
    if (alignment > 0.99 && minimumComponent > pathRayEpsilon) {
        pathAccumulateRadiance(
            sampleRadiance * (sampleThroughput / baseColor),
            sampleWeight,
            diffuseRadiance);
        return;
    }

    const vec3 perpendicularVector = cross(normalizedBase, white);
    if (dot(perpendicularVector, perpendicularVector) <=
        pathRayEpsilon * pathRayEpsilon) {
        pathAccumulateRadiance(
            sampleRadiance * sampleThroughput,
            sampleWeight,
            specularRadiance);
        return;
    }
    const vec3 perpendicular =
        safeNormalize(perpendicularVector, vec3(0.0));
    const vec3 planarThroughput =
        sampleThroughput -
        perpendicular * dot(perpendicular, sampleThroughput);
    const float sumDenominator = 1.0 + alignment;
    const float differenceDenominator = 1.0 - alignment;
    if (abs(sumDenominator) <= pathRayEpsilon ||
        abs(differenceDenominator) <= pathRayEpsilon) {
        pathAccumulateRadiance(
            sampleRadiance * sampleThroughput,
            sampleWeight,
            specularRadiance);
        return;
    }
    const float projectionOnWhite = dot(planarThroughput, white);
    const float projectionOnBase = dot(planarThroughput, normalizedBase);
    const float sum =
        (projectionOnWhite + projectionOnBase) / sumDenominator;
    const float difference =
        (projectionOnWhite - projectionOnBase) / differenceDenominator;
    const float baseCoefficient = 0.5 * (sum - difference);
    if (!isFiniteScalar(baseCoefficient)) {
        return;
    }
    const vec3 diffuseThroughput = baseCoefficient * normalizedBase;
    pathAccumulateRadiance(
        sampleRadiance * (baseCoefficient / baseLength),
        sampleWeight,
        diffuseRadiance);
    pathAccumulateRadiance(
        sampleRadiance * (sampleThroughput - diffuseThroughput),
        sampleWeight,
        specularRadiance);
}

float pathTransparentReflectionProbability(
    PathSurface surface, vec3 incoming, PathMediumStack media) {
    if (surface.opacity >= pathRayEpsilon) {
        return 1.0;
    }
    const float fresnel = pathPreciseFresnel(surface, incoming);
    const float reflectionProbability = media.count == 1U
        ? pathBaseReflectionProbability +
              (1.0 - pathBaseReflectionProbability) * fresnel
        : fresnel;
    return clamp(reflectionProbability, 0.0, 1.0);
}

void pathSampleDirectFromPrimary(
    PathSurface primarySurface,
    vec3 primaryRayDirection,
    float directProbability,
    float sampleWeight,
    uint pixelIndex,
    uint sampleIndex,
    inout vec4 directDiffuse,
    inout vec4 directSpecular) {
    if (!isFiniteVector(primarySurface.position) ||
        length(primarySurface.emission) > 0.0 ||
        !(directProbability > 0.0) || !(sampleWeight > 0.0)) {
        return;
    }
    PathSurface directSurface = primarySurface;
    if (directSurface.opacity < pathRayEpsilon) {
        directSurface.eta = 1.0 / max(directSurface.eta, 1.0e-4);
    }
    vec3 lightThroughput;
    vec3 lightRadiance;
    if (!pathSampleDirectLight(
            directSurface,
            -primaryRayDirection,
            pixelIndex,
            sampleIndex,
            0U,
            0U,
            0U,
            lightThroughput,
            lightRadiance)) {
        return;
    }
    lightThroughput /= max(directProbability, 1.0e-6);
    pathAccumulateInwardRadiance(
        primarySurface.baseColor,
        lightRadiance,
        lightThroughput,
        sampleWeight,
        directDiffuse,
        directSpecular);
}

void pathAccumulateIndirectEndpoint(
    PathSurface primarySurface,
    vec3 firstThroughput,
    vec3 downstreamRadiance,
    float sampleWeight,
    inout vec4 indirectDiffuse,
    inout vec4 indirectSpecular) {
    pathAccumulateInwardRadiance(
        primarySurface.opacity < pathRayEpsilon
            ? vec3(0.0)
            : primarySurface.baseColor,
        downstreamRadiance,
        firstThroughput,
        sampleWeight,
        indirectDiffuse,
        indirectSpecular);
}

void pathTraceIndirectFromPrimary(
    PathSurface primarySurface,
    vec3 primaryRayDirection,
    PathRayDifferential primaryDifferential,
    vec3 primaryPositionDx,
    vec3 primaryPositionDy,
    float indirectProbability,
    float sampleWeight,
    uint pixelIndex,
    uint sampleIndex,
    uint replicateIndex,
    inout vec4 indirectDiffuse,
    inout vec4 indirectSpecular,
    inout uint directSampleCount) {
    if (length(primarySurface.emission) > 0.0 ||
        !(indirectProbability > 0.0) || !(sampleWeight > 0.0)) {
        return;
    }

    PathMediumStack media;
    pathResetMediumStack(media);
    PathSurface firstSurface = primarySurface;
    const float firstAbsoluteIor = firstSurface.eta;
    if (firstSurface.opacity < pathRayEpsilon) {
        firstSurface.eta = 1.0 / max(firstSurface.eta, 1.0e-4);
    }
    const vec3 firstIncoming = -primaryRayDirection;
    const float firstReflectionProbability =
        pathTransparentReflectionProbability(firstSurface, firstIncoming, media);
    const float firstBranchRandom = pathRandom(
        pixelIndex,
        sampleIndex,
        0U,
        replicateIndex,
        pathDimensionBranch);
    const bool firstReflect = firstReflectionProbability >= 1.0 ||
        (firstReflectionProbability > 0.0 &&
         firstBranchRandom < firstReflectionProbability);
    vec3 rayDirection;
    vec3 firstThroughput;
    PathRayDifferential differential = primaryDifferential;
    bool validScatter = false;
    if (firstReflect) {
        validScatter = pathSampleReflection(
            firstSurface,
            firstIncoming,
            pathRandom(
                pixelIndex,
                sampleIndex,
                0U,
                replicateIndex,
                pathDimensionBsdfTechnique),
            pathRandom(
                pixelIndex,
                sampleIndex,
                0U,
                replicateIndex,
                pathDimensionBsdfFirst),
            pathRandom(
                pixelIndex,
                sampleIndex,
                0U,
                replicateIndex,
                pathDimensionBsdfSecond),
            rayDirection,
            firstThroughput);
        firstThroughput /= max(firstReflectionProbability, 1.0e-6);
        differential = pathReflectionDifferential(
            primaryRayDirection,
            firstSurface.surfaceNormal,
            primaryPositionDx,
            primaryPositionDy,
            primaryDifferential);
    } else {
        validScatter = pathSampleTransmission(
            firstSurface,
            firstIncoming,
            pathRandom(
                pixelIndex,
                sampleIndex,
                0U,
                replicateIndex,
                pathDimensionBsdfFirst),
            pathRandom(
                pixelIndex,
                sampleIndex,
                0U,
                replicateIndex,
                pathDimensionBsdfSecond),
            rayDirection,
            firstThroughput);
        firstThroughput /= max(1.0 - firstReflectionProbability, 1.0e-6);
        pathApplyTransmissionBoundary(media, firstSurface, firstAbsoluteIor);
    }
    firstThroughput /= max(indirectProbability, 1.0e-6);
    if (!validScatter || !isFiniteVector(rayDirection) ||
        !isFiniteVector(firstThroughput)) {
        return;
    }

    vec3 rayOrigin = primarySurface.position;
    vec3 downstreamThroughput = vec3(1.0);
    float roughnessFloor =
        primarySurface.roughness * pathRegularizationFactor;
    bool excludeDirectLight = true;

    for (uint depth = 1U; depth <= pathMaximumDepth; ++depth) {
        const PackedTraceHit hit = tracePackedScene(
            rayOrigin,
            pathParameters.cameraDirectionAndRayMinimum.w,
            rayDirection,
            rayMaximum,
            0U);
        if (!hit.hasHit) {
            if (!excludeDirectLight && pathHasEnvironment()) {
                pathAccumulateIndirectEndpoint(
                    primarySurface,
                    firstThroughput,
                    downstreamThroughput *
                        pathEnvironmentRadiance(rayDirection),
                    sampleWeight,
                    indirectDiffuse,
                    indirectSpecular);
            }
            return;
        }

        PathSurface surface;
        vec3 positionDx;
        vec3 positionDy;
        if (!pathPopulateSurface(
                hit,
                rayOrigin,
                rayDirection,
                differential,
                positionDx,
                positionDy,
                surface)) {
            return;
        }
        const float absoluteIor = surface.eta;
        roughnessFloor = max(
            roughnessFloor,
            pathRegularizationFactor * surface.roughness);
        surface.roughness = max(surface.roughness, roughnessFloor);
        const float rouletteProbability = clamp(
            1.0 - 0.5 * pathSafeSqrt(surface.roughness),
            0.5,
            1.0);
        bool maySampleDirect = rouletteProbability < 0.9;
        const vec3 absorption = pathMediumTransmission(
            pathCombinedAbsorption(media), hit.distance);

        if (length(surface.emission) > pathRayEpsilon) {
            if (!excludeDirectLight) {
                pathAccumulateIndirectEndpoint(
                    primarySurface,
                    firstThroughput,
                    downstreamThroughput * absorption * surface.emission,
                    sampleWeight,
                    indirectDiffuse,
                    indirectSpecular);
            }
            return;
        }

        if (surface.opacity < pathRayEpsilon) {
            surface.eta = pathRelativeEta(media, surface);
        }
        const vec3 incoming = -rayDirection;
        const float reflectionProbability =
            pathTransparentReflectionProbability(surface, incoming, media);
        const float branchRandom = pathRandom(
            pixelIndex,
            sampleIndex,
            depth,
            replicateIndex,
            pathDimensionBranch);
        const bool reflect = reflectionProbability >= 1.0 ||
            (reflectionProbability > 0.0 &&
             branchRandom < reflectionProbability);

        if (reflect) {
            maySampleDirect = maySampleDirect && surface.entering &&
                              media.count == 1U;
        } else {
            maySampleDirect = maySampleDirect && !surface.entering &&
                              media.count == 2U;
        }
        const bool terminateWithDirect = maySampleDirect &&
            (depth >= pathMaximumDepth ||
             pathRandom(
                 pixelIndex,
                 sampleIndex,
                 depth,
                 replicateIndex,
                 pathDimensionContinuation) > rouletteProbability);
        if (terminateWithDirect) {
            ++directSampleCount;
            const uint lightSampleCount = pathLightSampleCounts[depth];
            const float branchProbability = reflect
                                                ? reflectionProbability
                                                : 1.0 - reflectionProbability;
            float branchScale = 1.0 / max(branchProbability, 1.0e-6);
            if (depth < pathMaximumDepth) {
                branchScale /= max(1.0 - rouletteProbability, 1.0e-6);
            }
            for (uint lightIndex = 0U;
                 lightIndex < lightSampleCount;
                 ++lightIndex) {
                vec3 lightThroughput;
                vec3 lightRadiance;
                if (!pathSampleDirectLight(
                        surface,
                        incoming,
                        pixelIndex,
                        sampleIndex,
                        depth,
                        replicateIndex,
                        lightIndex,
                        lightThroughput,
                        lightRadiance)) {
                    continue;
                }
                vec3 downstreamRadiance =
                    downstreamThroughput * lightRadiance *
                    (lightThroughput * branchScale);
                if (!reflect) {
                    downstreamRadiance *= absorption;
                }
                pathAccumulateIndirectEndpoint(
                    primarySurface,
                    firstThroughput,
                    downstreamRadiance,
                    sampleWeight / float(lightSampleCount),
                    indirectDiffuse,
                    indirectSpecular);
            }
            return;
        }

        vec3 newDirection;
        vec3 scatterThroughput;
        bool scatterValid = false;
        if (reflect) {
            scatterValid = pathSampleReflection(
                surface,
                incoming,
                pathRandom(
                    pixelIndex,
                    sampleIndex,
                    depth,
                    replicateIndex,
                    pathDimensionBsdfTechnique),
                pathRandom(
                    pixelIndex,
                    sampleIndex,
                    depth,
                    replicateIndex,
                    pathDimensionBsdfFirst),
                pathRandom(
                    pixelIndex,
                    sampleIndex,
                    depth,
                    replicateIndex,
                    pathDimensionBsdfSecond),
                newDirection,
                scatterThroughput);
            scatterThroughput /= max(reflectionProbability, 1.0e-6);
            if (maySampleDirect) {
                scatterThroughput /= max(rouletteProbability, 1.0e-6);
            }
            differential = pathReflectionDifferential(
                rayDirection,
                surface.surfaceNormal,
                positionDx,
                positionDy,
                differential);
        } else {
            scatterValid = pathSampleTransmission(
                surface,
                incoming,
                pathRandom(
                    pixelIndex,
                    sampleIndex,
                    depth,
                    replicateIndex,
                    pathDimensionBsdfFirst),
                pathRandom(
                    pixelIndex,
                    sampleIndex,
                    depth,
                    replicateIndex,
                    pathDimensionBsdfSecond),
                newDirection,
                scatterThroughput);
            scatterThroughput /= max(1.0 - reflectionProbability, 1.0e-6);
            if (maySampleDirect) {
                scatterThroughput /= max(rouletteProbability, 1.0e-6);
            }
            pathApplyTransmissionBoundary(media, surface, absoluteIor);
        }
        if (!scatterValid || !isFiniteVector(newDirection) ||
            !isFiniteVector(scatterThroughput) ||
            depth >= pathMaximumDepth) {
            return;
        }
        downstreamThroughput *= scatterThroughput * absorption;
        if (!isFiniteVector(downstreamThroughput)) {
            return;
        }
        rayOrigin = surface.position;
        rayDirection = newDirection;
        excludeDirectLight = maySampleDirect;
    }
}

void raym0nadePathTraceInvocation() {
    const uvec2 tilePixel = gl_GlobalInvocationID.xy;
    const uvec2 tileExtent = pathParameters.tileExtentAndSampleRange.xy;
    if (tilePixel.x >= tileExtent.x || tilePixel.y >= tileExtent.y) {
        return;
    }
    const uvec2 imageExtent = pathParameters.imageExtentAndTileOrigin.xy;
    const uvec2 imagePixel =
        pathParameters.imageExtentAndTileOrigin.zw + tilePixel;
    if (imagePixel.x >= imageExtent.x || imagePixel.y >= imageExtent.y) {
        return;
    }
    const uint outputIndex = tilePixel.y * tileExtent.x + tilePixel.x;
    const uint pixelIndex = imagePixel.y * imageExtent.x + imagePixel.x;
    const uint sampleBase = pathParameters.tileExtentAndSampleRange.z;
    const uint sampleCount = pathParameters.tileExtentAndSampleRange.w;
    if (sampleBase == 0U) {
        pathResetOutput(outputIndex);
    }
    if (pathParameters.environmentFlagsTotalSppSceneVersionAndFlags.z !=
        pathPackedSceneVersion) {
        return;
    }

    const float imageX = float(imagePixel.x) - float(imageExtent.x) * 0.5;
    const float imageY = float(imagePixel.y) - float(imageExtent.y) * 0.5;
    const vec3 unnormalizedDirection =
        pathParameters.cameraDirectionAndRayMinimum.xyz +
        pathParameters.cameraPositionAndPixelScale.w *
            (imageX * pathParameters.cameraRightAndSampleWeight.xyz +
             imageY * pathParameters.cameraUpAndDirectProbability.xyz);
    const vec3 primaryRayDirection = normalizedOrZero(unnormalizedDirection);
    if (dot(primaryRayDirection, primaryRayDirection) == 0.0) {
        return;
    }
    const PathRayDifferential primaryDifferential =
        pathInitialDifferential(unnormalizedDirection);

    PathSurface primarySurface;
    vec3 primaryPositionDx = vec3(0.0);
    vec3 primaryPositionDy = vec3(0.0);
    if (sampleBase == 0U) {
        const PackedTraceHit primaryHit = tracePackedScene(
            pathParameters.cameraPositionAndPixelScale.xyz,
            pathParameters.cameraDirectionAndRayMinimum.w,
            primaryRayDirection,
            rayMaximum,
            0U);
        if (!primaryHit.hasHit) {
            pathOutputBuffer.pixels[outputIndex].emissionAndEta =
                vec4(pathEnvironmentRadiance(primaryRayDirection), 1.0);
            return;
        }
        if (!pathPopulateSurface(
                primaryHit,
                pathParameters.cameraPositionAndPixelScale.xyz,
                primaryRayDirection,
                primaryDifferential,
                primaryPositionDx,
                primaryPositionDy,
                primarySurface)) {
            return;
        }
        pathStoreSurface(outputIndex, primarySurface);
    } else {
        if (pathOutputBuffer.pixels[outputIndex]
                .materialEnteringHitAndDirectSampleCount.z == 0U) {
            return;
        }
        primarySurface = pathLoadSurface(outputIndex);
        const float hitDistance = length(
            primarySurface.position -
            pathParameters.cameraPositionAndPixelScale.xyz);
        pathPositionDifferentials(
            primaryRayDirection,
            hitDistance,
            primarySurface.shapeNormal,
            primaryDifferential,
            primaryPositionDx,
            primaryPositionDy);
    }

    vec4 directDiffuse = sampleBase == 0U
                             ? vec4(0.0)
                             : pathOutputBuffer.pixels[outputIndex].directDiffuse;
    vec4 directSpecular = sampleBase == 0U
                              ? vec4(0.0)
                              : pathOutputBuffer.pixels[outputIndex].directSpecular;
    vec4 indirectDiffuse = sampleBase == 0U
                               ? vec4(0.0)
                               : pathOutputBuffer.pixels[outputIndex].indirectDiffuse;
    vec4 indirectSpecular = sampleBase == 0U
                                ? vec4(0.0)
                                : pathOutputBuffer.pixels[outputIndex].indirectSpecular;
    uint directSampleCount = sampleBase == 0U
                                 ? 0U
                                 : pathOutputBuffer.pixels[outputIndex]
                                       .materialEnteringHitAndDirectSampleCount.w;
    const float directProbability = clamp(
        pathParameters.cameraUpAndDirectProbability.w, 0.0, 1.0);
    const float indirectProbability = 1.0 - directProbability;
    const uint totalSpp =
        pathParameters.environmentFlagsTotalSppSceneVersionAndFlags.y;
    float sampleWeight = pathParameters.cameraRightAndSampleWeight.w;
    if (!(sampleWeight > 0.0) || !isFiniteScalar(sampleWeight)) {
        sampleWeight = totalSpp > 0U ? 1.0 / float(totalSpp) : 0.0;
    }
    const uint replicateCount =
        primarySurface.opacity > 1.0 - pathRayEpsilon
            ? 1U
            : pathTransparentReplicates;

    for (uint batchSample = 0U; batchSample < sampleCount; ++batchSample) {
        const uint sampleIndex = sampleBase + batchSample;
        const float techniqueRandom = pathRandom(
            pixelIndex,
            sampleIndex,
            0U,
            0U,
            pathDimensionTechnique);
        const bool chooseDirect = directProbability >= 1.0 ||
            (directProbability > 0.0 &&
             techniqueRandom < directProbability);
        if (chooseDirect) {
            ++directSampleCount;
            pathSampleDirectFromPrimary(
                primarySurface,
                primaryRayDirection,
                directProbability,
                sampleWeight,
                pixelIndex,
                sampleIndex,
                directDiffuse,
                directSpecular);
            continue;
        }
        if (!(indirectProbability > 0.0)) {
            continue;
        }
        for (uint replicateIndex = 0U;
             replicateIndex < replicateCount;
             ++replicateIndex) {
            pathTraceIndirectFromPrimary(
                primarySurface,
                primaryRayDirection,
                primaryDifferential,
                primaryPositionDx,
                primaryPositionDy,
                indirectProbability,
                sampleWeight / float(replicateCount),
                pixelIndex,
                sampleIndex,
                replicateIndex,
                indirectDiffuse,
                indirectSpecular,
                directSampleCount);
        }
    }

    pathOutputBuffer.pixels[outputIndex].directDiffuse = directDiffuse;
    pathOutputBuffer.pixels[outputIndex].directSpecular = directSpecular;
    pathOutputBuffer.pixels[outputIndex].indirectDiffuse = indirectDiffuse;
    pathOutputBuffer.pixels[outputIndex].indirectSpecular = indirectSpecular;
    pathOutputBuffer.pixels[outputIndex]
        .materialEnteringHitAndDirectSampleCount.w = directSampleCount;
}

#endif
