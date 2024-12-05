#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <string>
#include <assimp/material.h>
#include <FreeImage.h>
#include "geometry.h"

const int TextureIdForDiffuseColor = 1;
const int TextureIdForSpecularColor = 2;
const int TextureIdForAmbientColor = 3;
const int TextureIdForEmission = 4;

const float gammaMap[256] = {0.00000000f,0.00000508f,0.00002333f,0.00005692f,0.00010719f,0.00017512f,0.00026154f,0.00036714f,0.00049250f,0.00063818f,0.00080466f,0.00099237f,0.00120174f,0.00143313f,0.00168692f,0.00196342f,0.00226295f,0.00258583f,0.00293232f,0.00330270f,0.00369724f,0.00411618f,0.00455975f,0.00502820f,0.00552175f,0.00604059f,0.00658496f,0.00715504f,0.00775103f,0.00837312f,0.00902149f,0.00969633f,0.01039780f,0.01112608f,0.01188133f,0.01266372f,0.01347340f,0.01431052f,0.01517524f,0.01606770f,0.01698805f,0.01793643f,0.01891298f,0.01991784f,0.02095113f,0.02201299f,0.02310356f,0.02422294f,0.02537128f,0.02654869f,0.02775528f,0.02899119f,0.03025652f,0.03155139f,0.03287592f,0.03423021f,0.03561437f,0.03702852f,0.03847275f,0.03994717f,0.04145189f,0.04298702f,0.04455263f,0.04614884f,0.04777576f,0.04943346f,0.05112205f,0.05284163f,0.05459228f,0.05637410f,0.05818718f,0.06003161f,0.06190748f,0.06381487f,0.06575388f,0.06772459f,0.06972709f,0.07176145f,0.07382777f,0.07592613f,0.07805659f,0.08021926f,0.08241421f,0.08464151f,0.08690125f,0.08919351f,0.09151836f,0.09387587f,0.09626612f,0.09868920f,0.10114516f,0.10363410f,0.10615607f,0.10871115f,0.11129942f,0.11392093f,0.11657579f,0.11926404f,0.12198573f,0.12474097f,0.12752980f,0.13035230f,0.13320853f,0.13609856f,0.13902247f,0.14198031f,0.14497215f,0.14799805f,0.15105806f,0.15415229f,0.15728074f,0.16044353f,0.16364069f,0.16687229f,0.17013839f,0.17343906f,0.17677434f,0.18014431f,0.18354902f,0.18698853f,0.19046290f,0.19397219f,0.19751646f,0.20109574f,0.20471014f,0.20835967f,0.21204442f,0.21576442f,0.21951973f,0.22331043f,0.22713655f,0.23099814f,0.23489527f,0.23882800f,0.24279638f,0.24680044f,0.25084025f,0.25491586f,0.25902736f,0.26317474f,0.26735809f,0.27157745f,0.27583286f,0.28012440f,0.28445208f,0.28881600f,0.29321617f,0.29765266f,0.30212551f,0.30663478f,0.31118053f,0.31576276f,0.32038158f,0.32503697f,0.32972905f,0.33445781f,0.33922336f,0.34402567f,0.34886485f,0.35374093f,0.35865393f,0.36360392f,0.36859095f,0.37361506f,0.37867627f,0.38377467f,0.38891026f,0.39408314f,0.39929333f,0.40454084f,0.40982577f,0.41514811f,0.42050794f,0.42590532f,0.43134022f,0.43681276f,0.44232297f,0.44787085f,0.45345649f,0.45907992f,0.46474114f,0.47044027f,0.47617728f,0.48195225f,0.48776522f,0.49361622f,0.49950528f,0.50543249f,0.51139784f,0.51740140f,0.52344316f,0.52952325f,0.53564173f,0.54179847f,0.54799360f,0.55422717f,0.56049925f,0.56680983f,0.57315898f,0.57954675f,0.58597308f,0.59243816f,0.59894192f,0.60548443f,0.61206573f,0.61868584f,0.62534481f,0.63204271f,0.63877958f,0.64555538f,0.65237021f,0.65922409f,0.66611707f,0.67304921f,0.68002045f,0.68703091f,0.69408065f,0.70116961f,0.70829791f,0.71546555f,0.72267258f,0.72991902f,0.73720491f,0.74453032f,0.75189519f,0.75929970f,0.76674372f,0.77422744f,0.78175080f,0.78931385f,0.79691666f,0.80455923f,0.81224161f,0.81996381f,0.82772595f,0.83552790f,0.84336984f,0.85125178f,0.85917372f,0.86713564f,0.87513769f,0.88317984f,0.89126217f,0.89938462f,0.90754735f,0.91575027f,0.92399347f,0.93227696f,0.94060081f,0.94896507f,0.95736969f,0.96581477f,0.97430032f,0.98282641f,0.99139297f,1.00000000f};

struct ImageData{
    int width, height, channels;
    std::vector<uint8_t> data;
    ImageData() : width(0), height(0), channels(0) {}

    [[nodiscard]] glm::vec3 get3(float u, float v) const;
    [[nodiscard]] glm::vec4 get4(float u, float v) const;
    [[nodiscard]] glm::vec4 get4_with_gamma(float u, float v) const;
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool hasTransparentPart() const;
};

std::string urlDecode(const std::string &src);

class Material {
public:

    ImageData texture[AI_TEXTURE_TYPE_MAX + 1];

    std::string name;
    float shininess;
    float opacity;
    vec3 diffuseColor;
    vec3 specularColor;
    vec3 ambientColor;
    vec3 emission;
    bool isNameEnabled;
    bool isShininessEnabled;
    bool isOpacityEnabled;
    bool isDiffuseColorEnabled;
    bool isSpecularColorEnabled;
    bool isAmbientColorEnabled;
    bool isEmissionEnabled;
    bool hasTransparentPart;

    Material();

    [[nodiscard]] const ImageData& getImage(int index) const;

    static bool loadImageFromDDS(ImageData &imageData, const std::string& filename);
    static bool loadImageFromPNG(ImageData &imageData, const std::string& filename);
    static bool loadImageFromJPG(ImageData &imageData, const std::string& filename);
    void loadImageFromFile(int index, const std::string& filename);

    // New method to load material properties
    void loadMaterialProperties(const aiMaterial* aiMat);
    [[nodiscard]] vec4 getDiffuseColor(float u, float v) const;
    [[nodiscard]] vec4 getNormal(float u, float v) const;
};

#endif // MATERIAL_H