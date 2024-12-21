#include <vector>
#include <cstring>
#include "image.h"

const float EDGE_THRESHOLD_MIN = 0.0312f;
const float EDGE_THRESHOLD_MAX = 0.125f;
const float SUBPIXEL_QUALITY = 0.75f;
const int QUALITY = 12; // 增加采样质量

void Image::FXAA() {

    const vec3 *data = pixelarray;
    vec3 *output = new vec3[width * height];

    auto getLuminance = [](const vec3& rgb) -> float {
        return dot(rgb, RGB_Weight);
    };

    for(int y = 0; y < height; y++) {
        for(int x = 0; x < width; x++) {
            float lumaM = getLuminance(data[y * width + x]);
            float lumaN = y > 0 ? getLuminance(data[(y-1) * width + x]) : lumaM;
            float lumaS = y < height-1 ? getLuminance(data[(y+1) * width + x]) : lumaM;
            float lumaE = x < width-1 ? getLuminance(data[y * width + x+1]) : lumaM;
            float lumaW = x > 0 ? getLuminance(data[y * width + x-1]) : lumaM;
            
            // 计算局部对比度
            float rangeMin = std::min({lumaN, lumaS, lumaE, lumaW});
            float rangeMax = std::max({lumaN, lumaS, lumaE, lumaW});
            float range = rangeMax - rangeMin;
            
            // 如果对比度太低则跳过
            if(range < std::max(EDGE_THRESHOLD_MIN, rangeMax * EDGE_THRESHOLD_MAX)) {
                output[y * width + x] = data[y * width + x];
                continue;
            }
            
            // 采样对角线像素
            float lumaNW = (y > 0 && x > 0) ? getLuminance(data[(y-1) * width + x-1]) : lumaM;
            float lumaNE = (y > 0 && x < width-1) ? getLuminance(data[(y-1) * width + x+1]) : lumaM;
            float lumaSW = (y < height-1 && x > 0) ? getLuminance(data[(y+1) * width + x-1]) : lumaM;
            float lumaSE = (y < height-1 && x < width-1) ? getLuminance(data[(y+1) * width + x+1]) : lumaM;
            
            // 计算边缘混合因子
            float edgeHorz = std::abs((lumaNW + lumaW + lumaSW) - (lumaNE + lumaE + lumaSE)) * (1.0f/3.0f);
            float edgeVert = std::abs((lumaNW + lumaN + lumaNE) - (lumaSW + lumaS + lumaSE)) * (1.0f/3.0f);
            
            bool isHorizontal = edgeHorz >= edgeVert;
            
            // 计算梯度和步长
            float stepLength = isHorizontal ? 1.0f/float(width) : 1.0f/float(height);
            float gradientStep = std::clamp(
                (isHorizontal ? edgeHorz : edgeVert) / range, 
                -2.0f, 
                2.0f
            );
            
            // 沿梯度方向采样
            vec2 uv(float(x)/float(width), float(y)/float(height));
            vec3 finalColor = data[y * width + x];
            float bestDelta = 0.0f;
            
            for(int i = 0; i < QUALITY; i++) {
                vec2 offset(
                    isHorizontal ? 0.0f : gradientStep * stepLength * float(i + 1),
                    isHorizontal ? gradientStep * stepLength * float(i + 1) : 0.0f
                );
                
                vec2 sampleUv = uv + offset;
                if(sampleUv.x < 0.0f || sampleUv.x > 1.0f || sampleUv.y < 0.0f || sampleUv.y > 1.0f) {
                    continue;
                }
                
                int sx = int(sampleUv.x * width);
                int sy = int(sampleUv.y * height);
                sx = std::clamp(sx, 0, width-1);
                sy = std::clamp(sy, 0, height-1);

                vec3 sampleColor = data[sy * width + sx];
                float sampleLuma = getLuminance(sampleColor);
                float delta = std::abs(sampleLuma - lumaM);
                
                if(delta > bestDelta) {
                    bestDelta = delta;
                    finalColor = sampleColor;
                }
            }
            
            // 子像素抗锯齿处理
            float subPixelOffset = std::min(
                (std::abs(lumaN + lumaS - 2.0f * lumaM) * 2.0f +
                 std::abs(lumaE + lumaW - 2.0f * lumaM)) * 0.25f,
                1.0f
            );
            
            // 混合最终颜色
            output[y * width + x] = mix(
                data[y * width + x],
                finalColor,
                subPixelOffset * SUBPIXEL_QUALITY
            );
        }
    }

    std::memcpy(pixelarray, output, width * height * sizeof(vec3));
    delete[] output;
}