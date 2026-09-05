#ifndef RAYM0NADE_PATH_CONTRACT_GLSL
#define RAYM0NADE_PATH_CONTRACT_GLSL

struct PackedAreaLight {
    vec4 centerAndPower;
    uvec4 triangleRangeAndReserved;
};

struct PackedAreaLightTriangle {
    uvec4 vertexIdsAndMaterialId;
    vec4 areaProbabilityCdfAndReserved;
};

struct PackedEnvironmentRow {
    vec4 probabilityCdfSolidAngleAndReserved;
};

struct PackedEnvironmentTexel {
    vec4 radianceAndConditionalCdf;
};

struct PathOutputPixel {
    vec4 shapeNormalAndRoughness;
    vec4 surfaceNormalAndSpecular;
    vec4 positionAndMetallic;
    vec4 baseColorAndOpacity;
    vec4 emissionAndEta;
    uvec4 materialEnteringHitAndDirectSampleCount;
    vec4 directDiffuse;
    vec4 directSpecular;
    vec4 indirectDiffuse;
    vec4 indirectSpecular;
};

layout(set = 0, binding = 5, std430) buffer PathOutputBuffer {
    PathOutputPixel pixels[];
} pathOutputBuffer;

layout(set = 0, binding = 9, std430) readonly buffer AreaLightBuffer {
    PackedAreaLight areaLights[];
} areaLightBuffer;

layout(set = 0, binding = 10, std430) readonly buffer AreaLightTriangleBuffer {
    PackedAreaLightTriangle areaLightTriangles[];
} areaLightTriangleBuffer;

layout(set = 0, binding = 11, std430) readonly buffer EnvironmentRowBuffer {
    PackedEnvironmentRow environmentRows[];
} environmentRowBuffer;

layout(set = 0, binding = 12, std430) readonly buffer EnvironmentTexelBuffer {
    PackedEnvironmentTexel environmentTexels[];
} environmentTexelBuffer;

layout(push_constant, std430) uniform PathTraceParameters {
    vec4 cameraPositionAndPixelScale;
    vec4 cameraDirectionAndRayMinimum;
    vec4 cameraUpAndDirectProbability;
    vec4 cameraRightAndSampleWeight;
    uvec4 imageExtentAndTileOrigin;
    uvec4 tileExtentAndSampleRange;
    uvec4 seedLightCountAndEnvironmentWidthHeight;
    uvec4 environmentFlagsTotalSppSceneVersionAndFlags;
} pathParameters;

bool packedPrimitiveRemapRequired() {
    return pathParameters.environmentFlagsTotalSppSceneVersionAndFlags.w != 0U;
}

struct PathRayDifferential {
    vec3 positionDx;
    vec3 positionDy;
    vec3 directionDx;
    vec3 directionDy;
};

struct PathSurface {
    vec3 shapeNormal;
    vec3 surfaceNormal;
    vec3 position;
    vec3 baseColor;
    vec3 emission;
    float specular;
    float roughness;
    float metallic;
    float opacity;
    float eta;
    uint materialId;
    bool entering;
};

struct PathMediumStack {
    uint materialIds[17];
    float iors[17];
    vec3 absorptions[17];
    uint count;
};

const uint pathPackedSceneVersion = 4U;
const uint pathInvalidSceneId = 0xffffffffU;
const uint pathEnvironmentHasImportance = 1U << 0U;
const uint pathMaterialHasSpecularTexture = 1U << 2U;
const uint pathMaterialHasEmissiveTexture = 1U << 3U;
const uint pathMaterialHasNormalTexture = 1U << 4U;

const uint pathMaximumDepth = 16U;
const uint pathMaximumMediumCount = 17U;
const uint pathTransparentReplicates = 16U;
const float pathRegularizationFactor = 1.0;
const float pathBaseReflectionProbability = 0.24;
const float pathMinimumPdf = 1.0e-8;
const float pathMinimumLightDistanceSquared = 1.0e-6;
const float pathThroughputLuminanceLimit = 64.0;
const float pathRayEpsilon = 1.0e-4;
const vec3 pathLuminanceWeights = vec3(0.3, 0.6, 0.1);

// Random dimensions are immutable ABI within this shader. Replicates use disjoint
// 1024-word regions, while each direct-light sample owns an eight-word subrange.
const uint pathReplicateDimensionStride = 1024U;
const uint pathDimensionTechnique = 0U;
const uint pathDimensionBranch = 1U;
const uint pathDimensionContinuation = 2U;
const uint pathDimensionBsdfTechnique = 3U;
const uint pathDimensionBsdfFirst = 4U;
const uint pathDimensionBsdfSecond = 5U;
const uint pathDimensionLightBase = 32U;
const uint pathDimensionLightStride = 8U;
const uint pathDimensionLightCategory = 0U;
const uint pathDimensionLightObject = 1U;
const uint pathDimensionLightTriangle = 2U;
const uint pathDimensionLightFirst = 3U;
const uint pathDimensionLightSecond = 4U;

bool pathFinite(vec4 value) {
    return !any(isnan(value)) && !any(isinf(value));
}

float pathSafeSqrt(float value) {
    return isFiniteScalar(value) ? sqrt(max(value, 0.0)) : 0.0;
}

float pathSquare(float value) {
    return value * value;
}

uint pathRandomDimensionAddress(uint replicateIndex, uint dimension) {
    return replicateIndex * pathReplicateDimensionStride + dimension;
}

uvec4 pathRandomBlockForDimension(
    uint pixelIndex,
    uint sampleIndex,
    uint bounceIndex,
    uint replicateIndex,
    uint dimension) {
    return raym0nadeCounterRandomBlockForDimension(
        pathParameters.seedLightCountAndEnvironmentWidthHeight.x,
        pixelIndex,
        sampleIndex,
        bounceIndex,
        pathRandomDimensionAddress(replicateIndex, dimension));
}

float pathRandomFromBlock(
    uvec4 block, uint replicateIndex, uint dimension) {
    return raym0nadeCounterRandomBlockOpen01(
        block, pathRandomDimensionAddress(replicateIndex, dimension));
}

float pathRandom(
    uint pixelIndex,
    uint sampleIndex,
    uint bounceIndex,
    uint replicateIndex,
    uint dimension) {
    const uvec4 block = pathRandomBlockForDimension(
        pixelIndex, sampleIndex, bounceIndex, replicateIndex, dimension);
    return pathRandomFromBlock(block, replicateIndex, dimension);
}

uint pathLightDimension(uint lightSampleIndex, uint component) {
    return pathDimensionLightBase +
           lightSampleIndex * pathDimensionLightStride + component;
}

void pathLimitThroughput(inout vec3 throughput) {
    const float luminance = dot(throughput, pathLuminanceWeights);
    if (isFiniteScalar(luminance) && luminance > pathThroughputLuminanceLimit) {
        throughput *= pathThroughputLuminanceLimit / luminance;
    }
}

void pathResetMediumStack(out PathMediumStack media) {
    // Inactive slots are never read, and entering a medium initializes its slot
    // before increasing the active count.
    media.materialIds[0] = pathInvalidSceneId;
    media.iors[0] = 1.0;
    media.absorptions[0] = vec3(1.0);
    media.count = 1U;
}

float pathCurrentIor(PathMediumStack media) {
    if (media.count == 0U) {
        return 1.0;
    }
    return media.iors[min(media.count, pathMaximumMediumCount) - 1U];
}

float pathIorAfterExit(PathMediumStack media, uint materialId) {
    for (int index = int(min(media.count, pathMaximumMediumCount)) - 1;
         index >= 0;
         --index) {
        if (media.materialIds[index] == materialId) {
            return index == 0 ? 1.0 : media.iors[index - 1];
        }
    }
    return 1.0;
}

vec3 pathCombinedAbsorption(PathMediumStack media) {
    vec3 absorption = vec3(1.0);
    const uint count = min(media.count, pathMaximumMediumCount);
    for (uint index = 0U; index < count; ++index) {
        absorption *= media.absorptions[index];
    }
    return absorption;
}

void pathEnterMedium(
    inout PathMediumStack media,
    uint materialId,
    float ior,
    vec3 absorption) {
    if (media.count >= pathMaximumMediumCount) {
        return;
    }
    media.materialIds[media.count] = materialId;
    media.iors[media.count] = max(ior, 1.0e-4);
    media.absorptions[media.count] = absorption;
    ++media.count;
}

void pathExitMedium(inout PathMediumStack media, uint materialId) {
    for (int index = int(min(media.count, pathMaximumMediumCount)) - 1;
         index >= 0;
         --index) {
        if (media.materialIds[index] != materialId) {
            continue;
        }
        for (uint moveIndex = uint(index); moveIndex + 1U < media.count; ++moveIndex) {
            media.materialIds[moveIndex] = media.materialIds[moveIndex + 1U];
            media.iors[moveIndex] = media.iors[moveIndex + 1U];
            media.absorptions[moveIndex] = media.absorptions[moveIndex + 1U];
        }
        media.count = max(media.count, 1U) - 1U;
        return;
    }
}

void pathApplyTransmissionBoundary(
    inout PathMediumStack media,
    PathSurface surface,
    float absoluteIor) {
    if (surface.entering) {
        pathEnterMedium(
            media, surface.materialId, absoluteIor, surface.baseColor);
    } else {
        pathExitMedium(media, surface.materialId);
    }
}

float pathRelativeEta(PathMediumStack media, PathSurface surface) {
    const float from = pathCurrentIor(media);
    const float to = surface.entering
                         ? surface.eta
                         : pathIorAfterExit(media, surface.materialId);
    return from / max(to, 1.0e-4);
}

vec3 pathMediumTransmission(vec3 absorption, float distance) {
    if (!isFiniteVector(absorption) || !isFiniteScalar(distance)) {
        return vec3(1.0);
    }
    const vec3 bounded = clamp(absorption, vec3(0.0), vec3(1.0));
    if (all(greaterThanEqual(bounded, vec3(1.0 - pathRayEpsilon)))) {
        return vec3(1.0);
    }
    const float exponent = min(32.0 * max(distance, 0.0), 3.402823466e+38);
    const vec3 result = pow(bounded, vec3(exponent));
    return isFiniteVector(result) ? result : vec3(1.0);
}

#endif
