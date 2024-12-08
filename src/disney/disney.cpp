/*
 * Disney BSDF Shading Model
 *
 * Original implementation by Joe Schutte (2018), as part of the Selas rendering engine.
 * Original repository: https://github.com/schuttejoe/Selas
 *
 * Modifications made to enhance the BSDF implementation:
 * - [列出你的修改，例如：“添加了对双向散射分布函数的改进”]
 *
 * Copyright (c) 2018 Joe Schutte
 * Modifications Copyright (c) 2024 Meng Chu and Yuchen Yang
 *
 * The original code is licensed under the MIT License. See the LICENSE file or
 * visit the original repository for details: https://github.com/schuttejoe/Selas
 */

#include "disney/disney.h"
#include "disney/utils.h"
#include "disney/Fresnel.h"

//=============================================================================================================================
static void CalculateLobePdfs(const BRDF &surface,
                              float &pSpecular, float &pDiffuse, float &pClearcoat, float &pSpecTrans) {
    float metallicBRDF = surface.metallic;
    float specularBSDF = (1.0f - surface.metallic) * surface.specTrans;
    float dielectricBRDF = (1.0f - surface.specTrans) * (1.0f - surface.metallic);

    float specularWeight = metallicBRDF + dielectricBRDF;
    float transmissionWeight = specularBSDF;
    float diffuseWeight = dielectricBRDF;
    float clearcoatWeight = 1.0f * ReLU(surface.clearcoat);

    float norm = 1.0f / (specularWeight + transmissionWeight + diffuseWeight + clearcoatWeight);

    pSpecular = specularWeight * norm;
    pSpecTrans = transmissionWeight * norm;
    pDiffuse = diffuseWeight * norm;
    pClearcoat = clearcoatWeight * norm;
}

//=============================================================================================================================
static float ThinTransmissionRoughness(float ior, float roughness) {
    return std::max(0.0f, (0.65f * ior - 0.35f) * roughness);
}

//=============================================================================================================================
static void CalculateAnisotropicParams(float roughness, float anisotropic, float &ax, float &ay) {
    float aspect = std::sqrt(1.0f - 0.9f * anisotropic);
    ax = std::max(0.001f, roughness * roughness / aspect);
    ay = std::max(0.001f, roughness * roughness * aspect);
}

//=============================================================================================================================
static glm::vec3 CalculateTint(glm::vec3 baseColor) {
    float luminance = glm::dot(glm::vec3(0.3f, 0.6f, 1.0f), baseColor);
    return (luminance > 0.0f) ? baseColor * (1.0f / luminance) : glm::vec3(1.0f);
}

//=============================================================================================================================
// -- "generalized" Trowbridge-Reitz curve ungeneralized with a hard-coded exponent of 1
static float GTR1(float absDotHL, float a) {
    if (a >= 1) {
        return M_1_PI;
    }

    float a2 = a * a;
    return (a2 - 1.0f) / (M_PI * std::log2(a2) * (1.0f + (a2 - 1.0f) * absDotHL * absDotHL));
}

//=============================================================================================================================

static float
EvaluateDisneyClearcoat(float clearcoat, float alpha, const glm::vec3 &wo, const glm::vec3 &wm, const glm::vec3 &wi,
                        float &fPdfW, float &rPdfW) {
    if (clearcoat <= 0.0f) {
        return 0.0f;
    }

    float absDotNH = std::abs(std::cosf(wm.y));
    float absDotNL = std::abs(std::cosf(wi.y));
    float absDotNV = std::abs(std::cosf(wo.y));
    float dotHL = glm::dot(wm, wi);

    float d = GTR1(absDotNH, lerp(0.1f, 0.001f, alpha));
    float f = Fresnel::Schlick(0.04f, dotHL);
    float gl = Bsdf::SeparableSmithGGXG1(wi, 0.25f);
    float gv = Bsdf::SeparableSmithGGXG1(wo, 0.25f);

    fPdfW = d / (4.0f * std::abs(glm::dot(wo, wm)));
    rPdfW = d / (4.0f * std::abs(glm::dot(wi, wm)));

    return 0.25f * clearcoat * d * f * gl * gv;
}

//=============================================================================================================================
static glm::vec3
EvaluateSheen(const BRDF &surface, const glm::vec3 &wo, const glm::vec3 &wm, const glm::vec3 &wi) {
    if (surface.sheen <= 0.0f) {
        return glm::vec3{0.0f};
    }

    float dotHL = std::abs(glm::dot(wm, wi));

    glm::vec3 tint = CalculateTint(surface.baseColor);
    return surface.sheen * lerp(glm::vec3(1.0f), tint, surface.sheenTint) * Fresnel::SchlickWeight(dotHL);
}

//=============================================================================================================================
static glm::vec3
DisneyFresnel(const BRDF &surface, const glm::vec3 &wo, const glm::vec3 &wm, const glm::vec3 &wi) {
    float dotHV = glm::dot(wm, wo);

    glm::vec3 tint = CalculateTint(surface.baseColor);

    glm::vec3 R0 =
            Fresnel::SchlickR0FromRelativeIOR(surface.relativeIOR) * lerp(glm::vec3(1.0f), tint, surface.specularTint);
    R0 = lerp(R0, surface.baseColor, surface.metallic);

    float dielectricFresnel = Fresnel::Dielectric(dotHV, 1.0f, surface.ior);
    glm::vec3 metallicFresnel = Fresnel::Schlick(R0, glm::dot(wi, wm));

    return lerp(glm::vec3(dielectricFresnel), metallicFresnel, surface.metallic);
}

//=============================================================================================================================
static glm::vec3
EvaluateDisneyBRDF(const BRDF &surface, const glm::vec3 &wo, const glm::vec3 &wm, const glm::vec3 &wi,
                   float &fPdf, float &rPdf) {
    fPdf = 0.0f;
    rPdf = 0.0f;

    float dotNL = CosTheta(wi);
    float dotNV = CosTheta(wo);
    if (dotNL <= 0.0f || dotNV <= 0.0f) {
        return glm::vec3{0.0f};
    }

    float ax, ay;
    CalculateAnisotropicParams(surface.roughness, surface.anisotropic, ax, ay);

    float d = Bsdf::GgxAnisotropicD(wm, ax, ay);
    float gl = Bsdf::SeparableSmithGGXG1(wi, wm, ax, ay);
    float gv = Bsdf::SeparableSmithGGXG1(wo, wm, ax, ay);

    glm::vec3 f = DisneyFresnel(surface, wo, wm, wi);

    Bsdf::GgxVndfAnisotropicPdf(wi, wm, wo, ax, ay, fPdf, rPdf);
    fPdf *= (1.0f / (4 * AbsDot(wo, wm)));
    rPdf *= (1.0f / (4 * AbsDot(wi, wm)));

    return d * gl * gv * f / (4.0f * dotNL * dotNV);
}

//=============================================================================================================================
static bool SampleDisneyBRDF(Generator &gen, const BRDF &surface, glm::vec3 v, BsdfSample &sample) {
    glm::vec3 wo = glm::normalize(MatrixMultiply(v, surface.worldToTangent));

    // -- Calculate Anisotropic params
    float ax, ay;
    CalculateAnisotropicParams(surface.roughness, surface.anisotropic, ax, ay);

    // -- Sample visible distribution of normals
    float r0 = gen();
    float r1 = gen();
    glm::vec3 wm = Bsdf::SampleGgxVndfAnisotropic(wo, ax, ay, r0, r1);

    // -- Reflect over wm
    glm::vec3 wi = glm::normalize(Reflect(wm, wo));
    if (CosTheta(wi) <= 0.0f) {
        sample.forwardPdfW = 0.0f;
        sample.reversePdfW = 0.0f;
        sample.reflectance = glm::vec3{0.0f};
        sample.wi = glm::vec3{0.0f};
        return false;
    }

    glm::vec3 F = DisneyFresnel(surface, wo, wm, wi);
    
    float G1v = Bsdf::SeparableSmithGGXG1(wo, wm, ax, ay);
    glm::vec3 specular = G1v * F;

    sample.flags = SurfaceEventFlags::eScatterEvent;
    sample.reflectance = specular;
    sample.wi = glm::normalize(MatrixMultiply(wi, glm::transpose(surface.worldToTangent)));
    Bsdf::GgxVndfAnisotropicPdf(wi, wm, wo, ax, ay, sample.forwardPdfW, sample.reversePdfW);

    sample.forwardPdfW *= (1.0f / (4 * AbsDot(wo, wm)));
    sample.reversePdfW *= (1.0f / (4 * AbsDot(wi, wm)));

    return true;
}

//=============================================================================================================================
static glm::vec3
EvaluateDisneySpecTransmission(const BRDF &surface, const glm::vec3 &wo, const glm::vec3 &wm,
                               const glm::vec3 &wi, float ax, float ay, bool thin) {
    float relativeIor = surface.relativeIOR;
    float n2 = relativeIor * relativeIor;

    float absDotNL = AbsCosTheta(wi);
    float absDotNV = AbsCosTheta(wo);
    float dotHL = Dot(wm, wi);
    float dotHV = Dot(wm, wo);
    float absDotHL = Absf(dotHL);
    float absDotHV = Absf(dotHV);

    float d = Bsdf::GgxAnisotropicD(wm, ax, ay);
    float gl = Bsdf::SeparableSmithGGXG1(wi, wm, ax, ay);
    float gv = Bsdf::SeparableSmithGGXG1(wo, wm, ax, ay);

    float f = Fresnel::Dielectric(dotHV, 1.0f, surface.ior);

    glm::vec3 color;
    if (thin)
        color = Sqrt(surface.baseColor);
    else
        color = surface.baseColor;

    // Note that we are intentionally leaving out the 1/n2 spreading factor since for VCM we will be evaluating particles with
    // this. That means we'll need to model the air-[other medium] transmission if we ever place the camera inside a non-air
    // medium.
    float c = (absDotHL * absDotHV) / (absDotNL * absDotNV);
    float t = (n2 / Square(dotHL + relativeIor * dotHV));
    return color * c * t * (1.0f - f) * gl * gv * d;
}

//=============================================================================================================================
static float EvaluateDisneyRetroDiffuse(const BRDF &surface, const glm::vec3 &wo, const glm::vec3 &wm,
                                        const glm::vec3 &wi) {
    float dotNL = AbsCosTheta(wi);
    float dotNV = AbsCosTheta(wo);

    float roughness = surface.roughness * surface.roughness;

    float rr = 0.5f + 2.0f * dotNL * dotNL * roughness;
    float fl = Fresnel::SchlickWeight(dotNL);
    float fv = Fresnel::SchlickWeight(dotNV);

    return rr * (fl + fv + fl * fv * (rr - 1.0f));
}

//=============================================================================================================================
static float
EvaluateDisneyDiffuse(const BRDF &surface, const glm::vec3 &wo, const glm::vec3 &wm, const glm::vec3 &wi,
                      bool thin) {
    float dotNL = AbsCosTheta(wi);
    float dotNV = AbsCosTheta(wo);

    float fl = Fresnel::SchlickWeight(dotNL);
    float fv = Fresnel::SchlickWeight(dotNV);

    float hanrahanKrueger = 0.0f;

    if (thin && surface.flatness > 0.0f) {
        float roughness = surface.roughness * surface.roughness;

        float dotHL = Dot(wm, wi);
        float fss90 = dotHL * dotHL * roughness;
        float fss = lerp(1.0f, fss90, fl) * lerp(1.0f, fss90, fv);

        float ss = 1.25f * (fss * (1.0f / (dotNL + dotNV) - 0.5f) + 0.5f);
        hanrahanKrueger = ss;
    }

    float lambert = 1.0f;
    float retro = EvaluateDisneyRetroDiffuse(surface, wo, wm, wi);
    float subsurfaceApprox = lerp(lambert, hanrahanKrueger, thin ? surface.flatness : 0.0f);

    return M_1_PI * (retro + subsurfaceApprox * (1.0f - 0.5f * fl) * (1.0f - 0.5f * fv));
}

//=============================================================================================================================
static bool
SampleDisneyClearcoat(Generator &gen, const BRDF &surface, const glm::vec3 &v, BsdfSample &sample) {
    glm::vec3 wo = MatrixMultiply(v, surface.worldToTangent);

    float a = 0.25f;
    float a2 = a * a;

    float r0 = gen();
    float r1 = gen();
    float cosTheta = std::sqrtf(std::max(0.0f, (1.0f - powf(a2, 1.0f - r0)) / (1.0f - a2)));
    float sinTheta = std::sqrtf(std::max(0.0f, 1.0f - cosTheta * cosTheta));
    float phi = 2.0 * M_PI * r1;

    glm::vec3 wm = glm::vec3(sinTheta * std::cosf(phi), cosTheta, sinTheta * std::sinf(phi));
    if (Dot(wm, wo) < 0.0f) {
        wm = -wm;
    }

    glm::vec3 wi = Reflect(wm, wo);
    if (Dot(wi, wo) < 0.0f) {
        return false;
    }

    float clearcoatWeight = surface.clearcoat;
    float clearcoatGloss = surface.clearcoatGloss;

    float dotNH = CosTheta(wm);
    float dotLH = Dot(wm, wi);

    float d = GTR1(Absf(dotNH), lerp(0.1f, 0.001f, clearcoatGloss));
    float f = Fresnel::Schlick(0.04f, dotLH);
    float g = Bsdf::SeparableSmithGGXG1(wi, 0.25f) * Bsdf::SeparableSmithGGXG1(wo, 0.25f);

    float fPdf = d / (4.0f * Dot(wo, wm));

    sample.reflectance = glm::vec3(0.25f * clearcoatWeight * g * f * d) / fPdf;
    sample.wi = glm::normalize(MatrixMultiply(wi, glm::transpose(surface.worldToTangent)));
    sample.forwardPdfW = fPdf;
    sample.reversePdfW = d / (4.0f * Dot(wi, wm));

    return true;
}

//=============================================================================================================================
static glm::vec3 CalculateExtinction(glm::vec3 apparantColor, float scatterDistance) {
    glm::vec3 a = apparantColor;
    glm::vec3 s = glm::vec3(1.9f) - a + 3.5f * (a - glm::vec3(0.8f)) * (a - glm::vec3(0.8f));

    return 1.0f / (s * scatterDistance);
}

//=============================================================================================================================
static bool SampleDisneySpecTransmission(Generator &gen, const BRDF &surface, glm::vec3 v, bool thin,
                                         BsdfSample &sample) {
    glm::vec3 wo = MatrixMultiply(v, surface.worldToTangent);
    if (CosTheta(wo) == 0.0) {
        sample.forwardPdfW = 0.0f;
        sample.reversePdfW = 0.0f;
        sample.reflectance = glm::vec3{0.0f};
        sample.wi = glm::vec3{0.0f};
        return false;
    }

    // -- Scale roughness based on IOR
    float rscaled = thin ? ThinTransmissionRoughness(surface.ior, surface.roughness) : surface.roughness;

    float tax, tay;
    CalculateAnisotropicParams(rscaled, surface.anisotropic, tax, tay);

    // -- Sample visible distribution of normals
    float r0 = gen();
    float r1 = gen();
    glm::vec3 wm = Bsdf::SampleGgxVndfAnisotropic(wo, tax, tay, r0, r1);

    float dotVH = Dot(wo, wm);
    if (wm.y < 0.0f) {
        dotVH = -dotVH;
    }

    float ni = wo.y > 0.0f ? 1.0f : surface.ior;
    float nt = wo.y > 0.0f ? surface.ior : 1.0f;
    float relativeIOR = ni / nt;

    // -- Disney uses the full dielectric Fresnel equation for transmission. We also importance sample F to switch between
    // -- refraction and reflection at glancing angles.
    float F = Fresnel::Dielectric(dotVH, 1.0f, surface.ior);

    // -- Since we're sampling the distribution of visible normals the pdf cancels out with a number of other terms.
    // -- We are left with the weight G2(wi, wo, wm) / G1(wi, wm) and since Disney uses a separable masking function
    // -- we get G1(wi, wm) * G1(wo, wm) / G1(wi, wm) = G1(wo, wm) as our weight.
    float G1v = Bsdf::SeparableSmithGGXG1(wo, wm, tax, tay);

    float pdf;

    glm::vec3 wi;
    if (gen() <= F) {
        wi = glm::normalize(Reflect(wm, wo));

        sample.flags = SurfaceEventFlags::eScatterEvent;
        sample.reflectance = G1v * surface.baseColor;

        float jacobian = (4 * AbsDot(wo, wm));
        pdf = F / jacobian;
    } else {
        if (thin) {
            // -- When the surface is thin so it refracts into and then out of the surface during this shading event.
            // -- So the ray is just reflected then flipped and we use the sqrt of the surface color.
            wi = Reflect(wm, wo);
            wi.y = -wi.y;
            sample.reflectance = G1v * Sqrt(surface.baseColor);

            // -- Since this is a thin surface we are not ending up inside of a volume so we treat this as a scatter event.
            sample.flags = SurfaceEventFlags::eScatterEvent;
        } else {
            if (Transmit(wm, wo, relativeIOR, wi)) {
                sample.flags = SurfaceEventFlags::eTransmissionEvent;
                sample.medium.phaseFunction =
                        dotVH > 0.0f ? MediumPhaseFunction::eIsotropic : MediumPhaseFunction::eVacuum;
                sample.medium.extinction = CalculateExtinction(surface.transmittanceColor, surface.scatterDistance);
            } else {
                sample.flags = SurfaceEventFlags::eScatterEvent;
                wi = Reflect(wm, wo);
            }

            sample.reflectance = G1v * surface.baseColor;
        }

        wi = glm::normalize(wi);

        float dotLH = Absf(Dot(wi, wm));
        float jacobian = dotLH / (Square(dotLH + surface.relativeIOR * dotVH));
        pdf = (1.0f - F) / jacobian;
    }

    if (CosTheta(wi) == 0.0f) {
        sample.forwardPdfW = 0.0f;
        sample.reversePdfW = 0.0f;
        sample.reflectance = glm::vec3{0.0f};
        sample.wi = glm::vec3{0.0f};
        return false;
    }

    if (surface.roughness < 0.01f) {
        sample.flags |= SurfaceEventFlags::eDiracEvent;
    }

    // -- calculate pdf terms
    Bsdf::GgxVndfAnisotropicPdf(wi, wm, wo, tax, tay, sample.forwardPdfW, sample.reversePdfW);
    sample.forwardPdfW *= pdf;
    sample.reversePdfW *= pdf;

    // -- convert wi back to world space
    sample.wi = glm::normalize(MatrixMultiply(wi, glm::transpose(surface.worldToTangent)));

    return true;
}

//=============================================================================================================================
static glm::vec3 SampleCosineWeightedHemisphere(float r0, float r1) {
    float r = std::sqrtf(r0);
    float theta = 2.0f * M_PI * r1;

    return glm::vec3(r * std::cosf(theta), std::sqrtf(std::fmax(0.0f, 1 - r0)), r * std::sinf(theta));
}

//=============================================================================================================================
static bool
SampleDisneyDiffuse(Generator &gen, const BRDF &surface, glm::vec3 v, bool thin, BsdfSample &sample) {
    glm::vec3 wo = glm::normalize(MatrixMultiply(v, surface.worldToTangent));

    float sign = Sign(CosTheta(wo));

    // -- Sample cosine lobe
    float r0 = gen();
    float r1 = gen();
    glm::vec3 wi = sign * SampleCosineWeightedHemisphere(r0, r1);
    glm::vec3 wm = glm::normalize(wi + wo);

    float dotNL = CosTheta(wi);
    if (dotNL == 0.0f) {
        sample.forwardPdfW = 0.0f;
        sample.reversePdfW = 0.0f;
        sample.reflectance = glm::vec3{0.0f};
        sample.wi = glm::vec3{0.0f};
        return false;
    }

    float dotNV = CosTheta(wo);

    float pdf;

    SurfaceEventFlags eventType = SurfaceEventFlags::eScatterEvent;

    glm::vec3 color = surface.baseColor;

    float p = gen();
    if (p <= surface.diffTrans) {
        wi = -wi;
        pdf = surface.diffTrans;

        if (thin)
            color = Sqrt(color);
        else {
            eventType = SurfaceEventFlags::eTransmissionEvent;

            sample.medium.phaseFunction = MediumPhaseFunction::eIsotropic;
            sample.medium.extinction = CalculateExtinction(surface.transmittanceColor, surface.scatterDistance);
        }
    } else {
        pdf = (1.0f - surface.diffTrans);
    }

    glm::vec3 sheen = EvaluateSheen(surface, wo, wm, wi);

    float diffuse = EvaluateDisneyDiffuse(surface, wo, wm, wi, thin);

    // Assert_(pdf > 0.0f);
    sample.reflectance = sheen + color * (diffuse / pdf);
    sample.wi = glm::normalize(MatrixMultiply(wi, glm::transpose(surface.worldToTangent)));
    sample.forwardPdfW = Absf(dotNL) * pdf;
    sample.reversePdfW = Absf(dotNV) * pdf;
    sample.flags = eventType;
    return true;
}

//=============================================================================================================================
glm::vec3 EvaluateDisney(const BRDF &surface, glm::vec3 v, glm::vec3 l, bool thin, float &forwardPdf,
                         float &reversePdf) {
    glm::vec3 wo = glm::normalize(MatrixMultiply(v, surface.worldToTangent));
    glm::vec3 wi = glm::normalize(MatrixMultiply(l, surface.worldToTangent));
    glm::vec3 wm = glm::normalize(wo + wi);

    float dotNV = CosTheta(wo);
    float dotNL = CosTheta(wi);

    glm::vec3 reflectance = glm::vec3{0.0f};
    forwardPdf = 0.0f;
    reversePdf = 0.0f;

    float pBRDF, pDiffuse, pClearcoat, pSpecTrans;
    CalculateLobePdfs(surface, pBRDF, pDiffuse, pClearcoat, pSpecTrans);

    float metallic = surface.metallic;
    float specTrans = surface.specTrans;

    // calculate all of the anisotropic params
    float ax, ay;
    CalculateAnisotropicParams(surface.roughness, surface.anisotropic, ax, ay);

    float diffuseWeight = (1.0f - metallic) * (1.0f - specTrans);
    float transWeight = (1.0f - metallic) * specTrans;

    // -- Clearcoat
    bool upperHemisphere = dotNL > 0.0f && dotNV > 0.0f;
    if (upperHemisphere && surface.clearcoat > 0.0f) {

        float forwardClearcoatPdfW;
        float reverseClearcoatPdfW;

        float clearcoat = EvaluateDisneyClearcoat(surface.clearcoat, surface.clearcoatGloss, wo, wm, wi,
                                                  forwardClearcoatPdfW, reverseClearcoatPdfW);
        reflectance += glm::vec3(clearcoat);
        forwardPdf += pClearcoat * forwardClearcoatPdfW;
        reversePdf += pClearcoat * reverseClearcoatPdfW;
    }

    // -- Diffuse
    if (diffuseWeight > 0.0f) {
        float forwardDiffusePdfW = AbsCosTheta(wi);
        float reverseDiffusePdfW = AbsCosTheta(wo);
        float diffuse = EvaluateDisneyDiffuse(surface, wo, wm, wi, thin);

        glm::vec3 sheen = EvaluateSheen(surface, wo, wm, wi);

        reflectance += diffuseWeight * (diffuse * surface.baseColor + sheen);

        forwardPdf += pDiffuse * forwardDiffusePdfW;
        reversePdf += pDiffuse * reverseDiffusePdfW;
    }

    // -- transmission
    if (transWeight > 0.0f) {

        // Scale roughness based on IOR (Burley 2015, Figure 15).
        float rscaled = thin ? ThinTransmissionRoughness(surface.ior, surface.roughness) : surface.roughness;
        float tax, tay;
        CalculateAnisotropicParams(rscaled, surface.anisotropic, tax, tay);

        glm::vec3 transmission = EvaluateDisneySpecTransmission(surface, wo, wm, wi, tax, tay, thin);
        reflectance += transWeight * transmission;

        float forwardTransmissivePdfW;
        float reverseTransmissivePdfW;
        Bsdf::GgxVndfAnisotropicPdf(wi, wm, wo, tax, tay, forwardTransmissivePdfW, reverseTransmissivePdfW);

        float dotLH = Dot(wm, wi);
        float dotVH = Dot(wm, wo);
        forwardPdf += pSpecTrans * forwardTransmissivePdfW / (Square(dotLH + surface.relativeIOR * dotVH));
        reversePdf += pSpecTrans * reverseTransmissivePdfW / (Square(dotVH + surface.relativeIOR * dotLH));
    }

    // -- specular
    if (upperHemisphere) {
        float forwardMetallicPdfW;
        float reverseMetallicPdfW;
        glm::vec3 specular = EvaluateDisneyBRDF(surface, wo, wm, wi, forwardMetallicPdfW, reverseMetallicPdfW);

        reflectance += specular;
        forwardPdf += pBRDF * forwardMetallicPdfW / (4 * AbsDot(wo, wm));
        reversePdf += pBRDF * reverseMetallicPdfW / (4 * AbsDot(wi, wm));
    }

    reflectance = reflectance * Absf(dotNL);

    return reflectance;
}

//=============================================================================================================================
bool SampleDisney(Generator &gen, const BRDF &surface, glm::vec3 v, bool thin, BsdfSample &sample) {
    float pSpecular;
    float pDiffuse;
    float pClearcoat;
    float pTransmission;
    CalculateLobePdfs(surface, pSpecular, pDiffuse, pClearcoat, pTransmission);

    bool success = false;

    float pLobe = 0.0f;
    float p = gen();
    if (p <= pSpecular) {
        success = SampleDisneyBRDF(gen, surface, v, sample);
        pLobe = pSpecular;
    } else if (p > pSpecular && p <= (pSpecular + pClearcoat)) {
        success = SampleDisneyClearcoat(gen, surface, v, sample);
        pLobe = pClearcoat;
    } else if (p > pSpecular + pClearcoat && p <= (pSpecular + pClearcoat + pDiffuse)) {
        success = SampleDisneyDiffuse(gen, surface, v, thin, sample);
        pLobe = pDiffuse;
    } else if (pTransmission >= 0.0f) {
        success = SampleDisneySpecTransmission(gen, surface, v, thin, sample);
        pLobe = pTransmission;
    } else {
        // -- Make sure we notice if this is occurring.
        sample.reflectance = glm::vec3(1000000.0f, 0.0f, 0.0f);
        sample.forwardPdfW = 0.000000001f;
        sample.reversePdfW = 0.000000001f;
    }

    if (pLobe > 0.0f) {
        sample.reflectance = sample.reflectance * (1.0f / pLobe);
        sample.forwardPdfW *= pLobe;
        sample.reversePdfW *= pLobe;
    }

    return success;
}