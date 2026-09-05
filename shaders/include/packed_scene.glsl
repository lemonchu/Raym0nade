#ifndef RAYM0NADE_PACKED_SCENE_GLSL
#define RAYM0NADE_PACKED_SCENE_GLSL

#extension GL_EXT_buffer_reference : require
#extension GL_EXT_buffer_reference2 : require

struct PackedVertex {
    vec4 positionAndNormalX;
    vec4 normalYZAndUv;
};

struct PackedMaterial {
    vec4 diffuseAndOpacity;
    vec4 emissionAndIor;
    vec4 transmissionAndRoughness;
    vec4 metallicSpecularAndReserved;
    uvec4 textureIds;
    uvec4 flagsAndReserved;
};

struct PackedTexture {
    uvec4 firstMipLevelCountAndExtent;
};

struct PackedTextureMip {
    uvec4 texelRangeAndExtent;
};

struct PackedTraceHit {
    bool hasHit;
    uint primitiveId;
    float distance;
    vec2 barycentrics;
};

layout(buffer_reference, std430, buffer_reference_align = 4) readonly buffer
    PackedTextureTexelPage {
    uint texelsRgba8[];
};

struct PackedTextureTexelPageRecord {
    PackedTextureTexelPage page;
    uint texelCount;
    uint reserved;
};

layout(set = 0, binding = 0) uniform accelerationStructureEXT topLevel;
layout(set = 0, binding = 1, std430) readonly buffer VertexBuffer {
    PackedVertex vertices[];
} vertexBuffer;
layout(set = 0, binding = 2, std430) readonly buffer IndexBuffer {
    uint indices[];
} indexBuffer;
layout(set = 0, binding = 3, std430) readonly buffer TriangleMaterialBuffer {
    uint materialIds[];
} triangleMaterialBuffer;
layout(set = 0, binding = 4, std430) readonly buffer MaterialBuffer {
    PackedMaterial materials[];
} materialBuffer;
layout(set = 0, binding = 6, std430) readonly buffer TextureBuffer {
    PackedTexture textures[];
} packedTextureBuffer;
layout(set = 0, binding = 7, std430) readonly buffer TextureMipBuffer {
    PackedTextureMip mipLevels[];
} textureMipBuffer;
layout(set = 0, binding = 8, std430) readonly buffer TextureTexelPageTableBuffer {
    // Page size in texels, page count, total texel count, and log2(page size).
    uvec4 paging;
    PackedTextureTexelPageRecord pages[];
} textureTexelPageTable;

const uint packedMaterialCutout = 1U << 0U;
const uint packedMaterialHasDiffuseTexture = 1U << 1U;
const float rayMaximum = 3.402823466e+38;
const float minimumNormalFloat = 1.175494351e-38;
const float cutoutAlphaThreshold = 1.0e-4;
const float diffuseGamma = 2.2;
const float pi = 3.14159265358979323846;

bool isFiniteScalar(float value) {
    return !isnan(value) && !isinf(value);
}

bool isFiniteVector(vec2 value) {
    return !any(isnan(value)) && !any(isinf(value));
}

bool isFiniteVector(vec3 value) {
    return !any(isnan(value)) && !any(isinf(value));
}

vec3 normalizedOrZero(vec3 value) {
    if (!isFiniteVector(value)) {
        return vec3(0.0);
    }
    const float scale = max(max(abs(value.x), abs(value.y)), abs(value.z));
    if (!(scale > minimumNormalFloat)) {
        return vec3(0.0);
    }
    const vec3 scaled = value / scale;
    const float lengthSquared = dot(scaled, scaled);
    if (!(lengthSquared > 0.0) || !isFiniteScalar(lengthSquared)) {
        return vec3(0.0);
    }
    return scaled * inversesqrt(lengthSquared);
}

vec3 safeNormalize(vec3 value, vec3 fallback) {
    const vec3 normalized = normalizedOrZero(value);
    if (dot(normalized, normalized) > 0.0) {
        return normalized;
    }
    return normalizedOrZero(fallback);
}

float wrapUnit(float value) {
    if (!isFiniteScalar(value)) {
        return 0.0;
    }
    const float wrapped = value - floor(value);
    return wrapped >= 0.0 && wrapped < 1.0 ? wrapped : 0.0;
}

vec4 unpackRgba8(uint value) {
    return vec4(
               float(value & 0xffU),
               float((value >> 8U) & 0xffU),
               float((value >> 16U) & 0xffU),
               float((value >> 24U) & 0xffU)) /
           255.0;
}

uint packedTextureTexel(uint texelIndex) {
    const uvec4 paging = textureTexelPageTable.paging;
    if (texelIndex >= paging.z || paging.x == 0U || paging.y == 0U ||
        paging.w >= 32U || (1U << paging.w) != paging.x) {
        return 0U;
    }
    const uint pageIndex = texelIndex >> paging.w;
    const uint pageOffset = texelIndex & (paging.x - 1U);
    if (pageIndex >= paging.y) {
        return 0U;
    }
    PackedTextureTexelPageRecord page =
        textureTexelPageTable.pages[pageIndex];
    if (pageOffset >= page.texelCount) {
        return 0U;
    }
    return page.page.texelsRgba8[pageOffset];
}

vec4 samplePackedMip(uint textureId, vec2 uv, uint relativeMipLevel) {
    if (!isFiniteVector(uv)) {
        return vec4(0.0);
    }
    const uvec4 texture =
        packedTextureBuffer.textures[textureId].firstMipLevelCountAndExtent;
    const uint mipLevel = min(relativeMipLevel, texture.y - 1U);
    const uvec4 mip =
        textureMipBuffer.mipLevels[texture.x + mipLevel].texelRangeAndExtent;
    const uvec2 extent = mip.zw;
    const vec2 wrappedUv = vec2(wrapUnit(uv.x), wrapUnit(1.0 - uv.y));
    const vec2 texelPosition = wrappedUv * vec2(extent);
    const vec2 flooredPosition = floor(texelPosition);
    const vec2 blend = texelPosition - flooredPosition;
    const uvec2 first = uvec2(flooredPosition) % extent;
    const uvec2 second = (first + uvec2(1U)) % extent;
    const uint rowStride = extent.x;
    const uint i00 = mip.x + first.y * rowStride + first.x;
    const uint i01 = mip.x + first.y * rowStride + second.x;
    const uint i10 = mip.x + second.y * rowStride + first.x;
    const uint i11 = mip.x + second.y * rowStride + second.x;
    const vec4 c00 = unpackRgba8(packedTextureTexel(i00));
    const vec4 c01 = unpackRgba8(packedTextureTexel(i01));
    const vec4 c10 = unpackRgba8(packedTextureTexel(i10));
    const vec4 c11 = unpackRgba8(packedTextureTexel(i11));
    const vec4 top = mix(c00, c01, blend.x);
    const vec4 bottom = mix(c10, c11, blend.x);
    return mix(top, bottom, blend.y);
}

vec4 samplePackedTexture(uint textureId, vec2 uv, float mipDepth) {
    const uvec4 texture =
        packedTextureBuffer.textures[textureId].firstMipLevelCountAndExtent;
    if (!isFiniteScalar(mipDepth)) {
        mipDepth = 0.0;
    }
    mipDepth = clamp(mipDepth, 0.0, float(texture.y - 1U));
    const uint firstLevel = uint(mipDepth);
    const uint secondLevel = min(firstLevel + 1U, texture.y - 1U);
    const float blend = mipDepth - float(firstLevel);
    return mix(
        samplePackedMip(textureId, uv, firstLevel),
        samplePackedMip(textureId, uv, secondLevel),
        blend);
}

uvec3 triangleVertexIds(uint primitiveId) {
    const uint indexOffset = primitiveId * 3U;
    return uvec3(
        indexBuffer.indices[indexOffset],
        indexBuffer.indices[indexOffset + 1U],
        indexBuffer.indices[indexOffset + 2U]);
}

vec3 vertexPosition(uint vertexId) {
    return vertexBuffer.vertices[vertexId].positionAndNormalX.xyz;
}

vec2 vertexUv(uint vertexId) {
    return vertexBuffer.vertices[vertexId].normalYZAndUv.zw;
}

vec3 barycentricWeights(vec2 rayQueryBarycentrics) {
    return vec3(
        1.0 - rayQueryBarycentrics.x - rayQueryBarycentrics.y,
        rayQueryBarycentrics.x,
        rayQueryBarycentrics.y);
}

vec2 interpolateUv(uint primitiveId, vec2 rayQueryBarycentrics) {
    const uvec3 vertexIds = triangleVertexIds(primitiveId);
    const vec3 weights = barycentricWeights(rayQueryBarycentrics);
    return weights.x * vertexUv(vertexIds.x) +
           weights.y * vertexUv(vertexIds.y) +
           weights.z * vertexUv(vertexIds.z);
}

bool candidateIsTransparentCutout(uint primitiveId, vec2 barycentrics) {
    const uint materialId = triangleMaterialBuffer.materialIds[primitiveId];
    const PackedMaterial material = materialBuffer.materials[materialId];
    if ((material.flagsAndReserved.x & packedMaterialCutout) == 0U) {
        return false;
    }
    const vec2 uv = interpolateUv(primitiveId, barycentrics);
    return samplePackedMip(material.textureIds.x, uv, 0U).a <
           cutoutAlphaThreshold;
}

PackedTraceHit tracePackedScene(
    vec3 origin, float tMinimum, vec3 direction, float tMaximum, uint rayFlags) {
    rayQueryEXT query;
    rayQueryInitializeEXT(
        query,
        topLevel,
        rayFlags,
        0xff,
        origin,
        tMinimum,
        direction,
        tMaximum);
    // Candidate order is implementation-defined, so a traversal-order cutout limit would
    // incorrectly turn ordinary rays into misses. Unlike the CPU's defensive 32-layer cap,
    // this query ignores every transparent candidate before accepting the nearest confirmed hit.
    while (rayQueryProceedEXT(query)) {
        if (rayQueryGetIntersectionTypeEXT(query, false) !=
            gl_RayQueryCandidateIntersectionTriangleEXT) {
            continue;
        }
        const uint primitiveId =
            rayQueryGetIntersectionPrimitiveIndexEXT(query, false);
        const vec2 barycentrics =
            rayQueryGetIntersectionBarycentricsEXT(query, false);
        if (!candidateIsTransparentCutout(primitiveId, barycentrics)) {
            rayQueryConfirmIntersectionEXT(query);
        }
    }
    if (rayQueryGetIntersectionTypeEXT(query, true) ==
        gl_RayQueryCommittedIntersectionNoneEXT) {
        return PackedTraceHit(false, 0U, 0.0, vec2(0.0));
    }
    return PackedTraceHit(
        true,
        rayQueryGetIntersectionPrimitiveIndexEXT(query, true),
        rayQueryGetIntersectionTEXT(query, true),
        rayQueryGetIntersectionBarycentricsEXT(query, true));
}

vec3 triangleShapeNormal(uint primitiveId, vec3 incomingDirection) {
    const uvec3 vertexIds = triangleVertexIds(primitiveId);
    const vec3 v0 = vertexPosition(vertexIds.x);
    const vec3 v1 = vertexPosition(vertexIds.y);
    const vec3 v2 = vertexPosition(vertexIds.z);
    vec3 shapeNormal =
        safeNormalize(cross(v1 - v0, v2 - v0), -incomingDirection);
    if (dot(shapeNormal, -incomingDirection) < 0.0) {
        shapeNormal = -shapeNormal;
    }
    return shapeNormal;
}

vec3 barycentricAtPoint(vec3 v0, vec3 v1, vec3 v2, vec3 point) {
    if (!isFiniteVector(v0) || !isFiniteVector(v1) ||
        !isFiniteVector(v2) || !isFiniteVector(point)) {
        return vec3(0.0);
    }
    const vec3 edge1 = v1 - v0;
    const vec3 edge2 = v2 - v0;
    const vec3 normal = cross(edge1, edge2);
    const float denominator = dot(normal, normal);
    if (!(denominator > minimumNormalFloat) ||
        !isFiniteScalar(denominator)) {
        return vec3(0.0);
    }
    const vec3 offset = point - v0;
    const float second = dot(cross(offset, edge2), normal) / denominator;
    const float third = dot(cross(edge1, offset), normal) / denominator;
    const vec3 result = vec3(1.0 - second - third, second, third);
    return isFiniteVector(result) ? result : vec3(0.0);
}

vec2 textureDerivative(uint primitiveId, vec3 positionDelta) {
    const uvec3 vertexIds = triangleVertexIds(primitiveId);
    const vec3 v0 = vertexPosition(vertexIds.x);
    const vec3 coordinates =
        barycentricAtPoint(
            v0,
            vertexPosition(vertexIds.y),
            vertexPosition(vertexIds.z),
            v0 + positionDelta);
    return (coordinates.x - 1.0) * vertexUv(vertexIds.x) +
           coordinates.y * vertexUv(vertexIds.y) +
           coordinates.z * vertexUv(vertexIds.z);
}

float diffuseMipDepth(uint textureId, float footprint) {
    const uint textureWidth =
        packedTextureBuffer.textures[textureId].firstMipLevelCountAndExtent.z;
    if (!isFiniteScalar(footprint) || !(footprint > 0.0) ||
        textureWidth == 0U) {
        return 0.0;
    }
    return max(0.0, log2(footprint * float(textureWidth)));
}

vec3 evaluateBaseColor(
    uint primitiveId, vec2 barycentrics, float textureFootprint) {
    const uint materialId = triangleMaterialBuffer.materialIds[primitiveId];
    const PackedMaterial material = materialBuffer.materials[materialId];
    if (material.diffuseAndOpacity.w < cutoutAlphaThreshold) {
        return material.transmissionAndRoughness.xyz;
    }
    if ((material.flagsAndReserved.x &
         packedMaterialHasDiffuseTexture) == 0U) {
        return max(material.diffuseAndOpacity.xyz, vec3(0.0));
    }
    const uint textureId = material.textureIds.x;
    const vec2 uv = interpolateUv(primitiveId, barycentrics);
    const float mipDepth =
        diffuseMipDepth(textureId, textureFootprint);
    const vec3 encoded =
        samplePackedTexture(textureId, uv, mipDepth).rgb;
    const vec3 linear = pow(max(encoded, vec3(0.0)), vec3(diffuseGamma));
    return linear * material.diffuseAndOpacity.xyz;
}

#endif
