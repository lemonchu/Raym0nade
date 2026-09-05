#ifndef RAYM0NADE_PATH_BSDF_GLSL
#define RAYM0NADE_PATH_BSDF_GLSL

float pathFifthPower(float value) {
    const float squared = value * value;
    return squared * squared * value;
}

float pathSchlickWeight(float cosine) {
    return pathFifthPower(clamp(1.0 - cosine, 0.0, 1.0));
}

float pathDielectricFresnel(float incidentCosine, float eta) {
    if (!isFiniteScalar(incidentCosine) || !isFiniteScalar(eta) || eta <= 0.0) {
        return 1.0;
    }
    incidentCosine = clamp(abs(incidentCosine), 0.0, 1.0);
    const float transmittedSineSquared =
        pathSquare(eta) * max(0.0, 1.0 - pathSquare(incidentCosine));
    if (transmittedSineSquared >= 1.0) {
        return 1.0;
    }
    const float transmittedCosine = pathSafeSqrt(1.0 - transmittedSineSquared);
    const float perpendicularDenominator =
        eta * incidentCosine + transmittedCosine;
    const float parallelDenominator =
        incidentCosine + eta * transmittedCosine;
    if (perpendicularDenominator <= pathMinimumPdf ||
        parallelDenominator <= pathMinimumPdf) {
        return 1.0;
    }
    const float perpendicular =
        (eta * incidentCosine - transmittedCosine) /
        perpendicularDenominator;
    const float parallel =
        (incidentCosine - eta * transmittedCosine) /
        parallelDenominator;
    return clamp(
        0.5 * (pathSquare(perpendicular) + pathSquare(parallel)),
        0.0,
        1.0);
}

float pathGtr1(float normalDotHalf, float alpha) {
    alpha = clamp(alpha, 1.0e-3, 0.9999);
    const float alphaSquared = alpha * alpha;
    const float denominator =
        pi * log(alphaSquared) *
        (1.0 + (alphaSquared - 1.0) * normalDotHalf * normalDotHalf);
    if (abs(denominator) <= pathMinimumPdf || !isFiniteScalar(denominator)) {
        return 0.0;
    }
    const float value = (alphaSquared - 1.0) / denominator;
    return isFiniteScalar(value) ? value : 0.0;
}

float pathGtr2(float normalDotHalf, float alpha) {
    alpha = max(alpha, 1.0e-3);
    const float alphaSquared = alpha * alpha;
    const float term =
        1.0 + (alphaSquared - 1.0) * normalDotHalf * normalDotHalf;
    const float value = alphaSquared /
                        max(pi * term * term, pathMinimumPdf);
    return isFiniteScalar(value) ? value : 0.0;
}

float pathSmithGgx(float normalDotDirection, float alpha) {
    if (normalDotDirection <= 0.0) {
        return 0.0;
    }
    const float alphaSquared = alpha * alpha;
    const float cosineSquared = normalDotDirection * normalDotDirection;
    const float root = pathSafeSqrt(
        alphaSquared + cosineSquared - alphaSquared * cosineSquared);
    return 1.0 / max(normalDotDirection + root, pathMinimumPdf);
}

float pathCosinePdf(vec3 normal, vec3 direction) {
    return max(0.0, dot(normal, direction)) / pi;
}

float pathMicrofacetReflectionPdf(
    vec3 normal, vec3 incoming, vec3 outgoing, float roughness) {
    const vec3 halfVector = safeNormalize(incoming + outgoing, vec3(0.0));
    if (dot(halfVector, halfVector) == 0.0) {
        return 0.0;
    }
    const float normalDotHalf = max(0.0, dot(normal, halfVector));
    const float outgoingDotHalf = abs(dot(outgoing, halfVector));
    if (normalDotHalf <= 0.0 || outgoingDotHalf <= pathMinimumPdf) {
        return 0.0;
    }
    const float pdf = pathGtr2(normalDotHalf, roughness) * normalDotHalf /
                      (4.0 * outgoingDotHalf);
    return isFiniteScalar(pdf) ? pdf : 0.0;
}

vec3 pathEvaluateReflection(
    PathSurface surface, vec3 incoming, vec3 outgoing) {
    const vec3 normal = surface.surfaceNormal;
    const float normalDotOutgoing = dot(normal, outgoing);
    const float normalDotIncoming = dot(normal, incoming);
    if (normalDotOutgoing <= 0.0 || normalDotIncoming <= 0.0) {
        return vec3(0.0);
    }
    const vec3 halfVector = safeNormalize(outgoing + incoming, vec3(0.0));
    if (dot(halfVector, halfVector) == 0.0) {
        return vec3(0.0);
    }
    const float normalDotHalf = max(0.0, dot(normal, halfVector));
    const float outgoingDotHalf = clamp(dot(outgoing, halfVector), 0.0, 1.0);
    const vec3 baseColor = max(surface.baseColor, vec3(0.0));
    const vec3 specularColor = mix(vec3(surface.specular), baseColor, surface.metallic);
    const float outgoingFresnel = pathSchlickWeight(normalDotOutgoing);
    const float incomingFresnel = pathSchlickWeight(normalDotIncoming);
    const float diffuseAtGrazing =
        0.5 + 2.0 * outgoingDotHalf * outgoingDotHalf * surface.roughness;
    const float diffuse =
        mix(1.0, diffuseAtGrazing, outgoingFresnel) *
        mix(1.0, diffuseAtGrazing, incomingFresnel);
    const bool transmissive = surface.opacity < pathRayEpsilon;
    const float fresnel = transmissive
                              ? pathDielectricFresnel(outgoingDotHalf, surface.eta)
                              : pathSchlickWeight(outgoingDotHalf);
    const vec3 specularFresnel = transmissive
                                     ? vec3(fresnel)
                                     : mix(specularColor, vec3(1.0), fresnel);
    const float distribution =
        pathGtr2(normalDotHalf, surface.roughness);
    const float geometry =
        pathSmithGgx(normalDotOutgoing, surface.roughness) *
        pathSmithGgx(normalDotIncoming, surface.roughness);
    const float clearcoatDistribution =
        pathGtr1(normalDotHalf, mix(0.1, 0.001, 0.2));
    const float clearcoatGeometry =
        pathSmithGgx(normalDotOutgoing, 0.25) *
        pathSmithGgx(normalDotIncoming, 0.25);
    const float clearcoatFresnel = mix(0.04, 1.0, fresnel);

    vec3 value = (1.0 / pi) * diffuse * baseColor *
                 (1.0 - surface.metallic) * surface.opacity;
    value += specularFresnel * distribution * geometry;
    value += surface.opacity * 0.375 * clearcoatGeometry *
             clearcoatFresnel * clearcoatDistribution * vec3(1.0);
    value *= normalDotOutgoing;
    return isFiniteVector(value) ? value : vec3(0.0);
}

vec3 pathEvaluateTransmission(
    PathSurface surface, vec3 incoming, vec3 outgoing) {
    const vec3 normal = surface.surfaceNormal;
    const float normalDotOutgoing = dot(normal, outgoing);
    const float normalDotIncoming = dot(normal, incoming);
    if (normalDotOutgoing >= 0.0 || normalDotIncoming <= 0.0) {
        return vec3(0.0);
    }
    vec3 halfVector = safeNormalize(outgoing + surface.eta * incoming, vec3(0.0));
    if (dot(halfVector, halfVector) == 0.0) {
        return vec3(0.0);
    }
    if (dot(normal, halfVector) < 0.0) {
        halfVector = -halfVector;
    }
    const float normalDotHalf = dot(normal, halfVector);
    const float outgoingDotHalf = dot(outgoing, halfVector);
    const float incomingDotHalf = dot(incoming, halfVector);
    if (normalDotHalf <= 0.0 ||
        outgoingDotHalf * normalDotOutgoing < 0.0 ||
        incomingDotHalf * normalDotIncoming < 0.0) {
        return vec3(0.0);
    }
    const float denominator =
        surface.eta * incomingDotHalf + outgoingDotHalf;
    if (abs(denominator) <= pathMinimumPdf || !isFiniteScalar(denominator)) {
        return vec3(0.0);
    }
    const float distribution =
        pathGtr2(normalDotHalf, surface.roughness);
    const float smithProduct =
        pathSmithGgx(abs(normalDotOutgoing), surface.roughness) *
        pathSmithGgx(abs(normalDotIncoming), surface.roughness);
    const float maskingShadowing =
        4.0 * abs(normalDotOutgoing * normalDotIncoming) * smithProduct;
    const float fresnel =
        pathDielectricFresnel(incomingDotHalf, surface.eta);
    float value = distribution * (1.0 - fresnel) * maskingShadowing *
                  pathSquare(surface.eta);
    value *= abs(outgoingDotHalf * incomingDotHalf) /
             max(abs(normalDotOutgoing * normalDotIncoming), pathMinimumPdf);
    value *= -normalDotOutgoing / (denominator * denominator);
    return isFiniteScalar(value) ? vec3(value) : vec3(0.0);
}

vec3 pathEvaluateBsdf(
    PathSurface surface, vec3 incoming, vec3 outgoing) {
    if (surface.opacity > 1.0 - pathRayEpsilon || surface.entering) {
        return pathEvaluateReflection(surface, incoming, outgoing);
    }
    return pathEvaluateTransmission(surface, incoming, outgoing);
}

vec3 pathRefractDirection(vec3 incoming, vec3 normal, float eta) {
    const float normalDotIncoming = dot(normal, incoming);
    const float discriminant =
        1.0 - eta * eta *
                  (1.0 - normalDotIncoming * normalDotIncoming);
    if (discriminant < 0.0 || !isFiniteScalar(discriminant)) {
        return vec3(uintBitsToFloat(0x7fc00000U));
    }
    return (eta * normalDotIncoming - pathSafeSqrt(discriminant)) * normal -
           eta * incoming;
}

float pathPreciseFresnel(PathSurface surface, vec3 incoming) {
    const float normalDotIncoming = dot(surface.surfaceNormal, incoming);
    if (normalDotIncoming < 0.0) {
        return 1.0;
    }
    const vec3 outgoing = pathRefractDirection(
        incoming, surface.surfaceNormal, surface.eta);
    if (!isFiniteVector(outgoing)) {
        return 1.0;
    }
    return pathDielectricFresnel(normalDotIncoming, surface.eta);
}

bool pathSampleReflection(
    PathSurface surface,
    vec3 incoming,
    float techniqueRandom,
    float firstRandom,
    float secondRandom,
    out vec3 outgoing,
    out vec3 throughput) {
    const bool sampleCosine = techniqueRandom < 0.5;
    vec3 tangent;
    vec3 bitangent;
    pathMakeTangentSpace(surface.surfaceNormal, incoming, tangent, bitangent);
    float selectedPdf = 0.0;
    if (sampleCosine) {
        const float radius = pathSafeSqrt(firstRandom);
        const float azimuth = secondRandom * 2.0 * pi;
        const float z = pathSafeSqrt(1.0 - radius * radius);
        outgoing = radius * cos(azimuth) * tangent +
                   radius * sin(azimuth) * bitangent +
                   z * surface.surfaceNormal;
        selectedPdf = pathCosinePdf(surface.surfaceNormal, outgoing);
    } else {
        const float azimuth = secondRandom * 2.0 * pi;
        const float roughnessSquared =
            pathSquare(max(surface.roughness, 1.0e-3));
        const float cosine = pathSafeSqrt(
            (1.0 - firstRandom) /
            (1.0 + (roughnessSquared - 1.0) * firstRandom));
        const float sine = pathSafeSqrt(1.0 - cosine * cosine);
        const vec3 halfVector =
            sine * cos(azimuth) * tangent +
            sine * sin(azimuth) * bitangent +
            cosine * surface.surfaceNormal;
        outgoing = 2.0 * dot(incoming, halfVector) * halfVector - incoming;
        selectedPdf = pathMicrofacetReflectionPdf(
            surface.surfaceNormal, incoming, outgoing, surface.roughness);
    }

    if (!isFiniteVector(outgoing) || selectedPdf <= 0.0 ||
        dot(outgoing, surface.shapeNormal) <= 0.0) {
        outgoing = vec3(0.0);
        throughput = vec3(0.0);
        return false;
    }
    const float diffusePdf = pathCosinePdf(surface.surfaceNormal, outgoing);
    const float specularPdf = pathMicrofacetReflectionPdf(
        surface.surfaceNormal, incoming, outgoing, surface.roughness);
    const float mixturePdf = 0.5 * (diffusePdf + specularPdf);
    if (mixturePdf <= pathMinimumPdf || !isFiniteScalar(mixturePdf)) {
        outgoing = vec3(0.0);
        throughput = vec3(0.0);
        return false;
    }
    throughput = pathEvaluateReflection(surface, incoming, outgoing) /
                 mixturePdf;
    pathLimitThroughput(throughput);
    return isFiniteVector(throughput);
}

bool pathSampleTransmission(
    PathSurface surface,
    vec3 incoming,
    float firstRandom,
    float secondRandom,
    out vec3 outgoing,
    out vec3 throughput) {
    if (abs(surface.eta - 1.0) < pathRayEpsilon) {
        outgoing = -incoming;
        throughput = vec3(1.0);
        return true;
    }
    vec3 tangent;
    vec3 bitangent;
    pathMakeTangentSpace(surface.surfaceNormal, incoming, tangent, bitangent);
    const float roughnessSquared =
        pathSquare(max(surface.roughness, 1.0e-3));
    const float cosine = pathSafeSqrt(
        (1.0 - firstRandom) /
        (1.0 + (roughnessSquared - 1.0) * firstRandom));
    const float sine = pathSafeSqrt(1.0 - cosine * cosine);
    const float azimuth = secondRandom * 2.0 * pi;
    const vec3 halfVector =
        sine * cos(azimuth) * tangent +
        sine * sin(azimuth) * bitangent +
        cosine * surface.surfaceNormal;
    outgoing = pathRefractDirection(incoming, halfVector, surface.eta);
    if (!isFiniteVector(outgoing) ||
        dot(outgoing, surface.shapeNormal) >= 0.0) {
        outgoing = vec3(0.0);
        throughput = vec3(0.0);
        return false;
    }
    const float normalDotHalf =
        max(0.0, dot(surface.surfaceNormal, halfVector));
    const float incomingDotHalf = dot(incoming, halfVector);
    const float outgoingDotHalf = dot(outgoing, halfVector);
    const float denominator =
        surface.eta * incomingDotHalf + outgoingDotHalf;
    if (abs(denominator) <= pathMinimumPdf || !isFiniteScalar(denominator)) {
        outgoing = vec3(0.0);
        throughput = vec3(0.0);
        return false;
    }
    const float halfPdf =
        pathGtr2(normalDotHalf, surface.roughness) * normalDotHalf;
    const float jacobian =
        abs(outgoingDotHalf) / (denominator * denominator);
    const float pdf = halfPdf * jacobian;
    if (pdf <= pathMinimumPdf || !isFiniteScalar(pdf)) {
        outgoing = vec3(0.0);
        throughput = vec3(0.0);
        return false;
    }
    throughput = pathEvaluateTransmission(surface, incoming, outgoing) / pdf;
    pathLimitThroughput(throughput);
    return isFiniteVector(throughput);
}

#endif
