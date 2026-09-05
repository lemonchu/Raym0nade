#ifndef RAYM0NADE_COUNTER_RNG_GLSL
#define RAYM0NADE_COUNTER_RNG_GLSL

const uint raym0nadePhiloxMultiplier0 = 0xD2511F53U;
const uint raym0nadePhiloxMultiplier1 = 0xCD9E8D57U;
const uint raym0nadePhiloxWeyl0 = 0x9E3779B9U;
const uint raym0nadePhiloxWeyl1 = 0xBB67AE85U;

uvec4 raym0nadePhiloxRound(uvec4 counter, uvec2 key) {
    uint firstHigh;
    uint firstLow;
    uint secondHigh;
    uint secondLow;
    umulExtended(raym0nadePhiloxMultiplier0, counter.x, firstHigh, firstLow);
    umulExtended(raym0nadePhiloxMultiplier1, counter.z, secondHigh, secondLow);
    return uvec4(
        secondHigh ^ counter.y ^ key.x,
        secondLow,
        firstHigh ^ counter.w ^ key.y,
        firstLow);
}

uvec4 raym0nadePhilox4x32TenRounds(uvec4 counter, uvec2 key) {
    for (int roundIndex = 0; roundIndex < 10; ++roundIndex) {
        counter = raym0nadePhiloxRound(counter, key);
        key += uvec2(raym0nadePhiloxWeyl0, raym0nadePhiloxWeyl1);
    }
    return counter;
}

uvec4 raym0nadeCounterRandomBlock(
    uint seed,
    uint pixelIndex,
    uint sampleIndex,
    uint bounceIndex,
    uint dimensionBlock) {
    return raym0nadePhilox4x32TenRounds(
        uvec4(pixelIndex, sampleIndex, bounceIndex, dimensionBlock),
        uvec2(seed, 0U));
}

uvec4 raym0nadeCounterRandomBlockForDimension(
    uint seed,
    uint pixelIndex,
    uint sampleIndex,
    uint bounceIndex,
    uint dimension) {
    return raym0nadeCounterRandomBlock(
        seed, pixelIndex, sampleIndex, bounceIndex, dimension >> 2U);
}

uint raym0nadeCounterRandomBlockWord(uvec4 block, uint dimension) {
    return block[dimension & 3U];
}

uint raym0nadeCounterRandomUint32(
    uint seed,
    uint pixelIndex,
    uint sampleIndex,
    uint bounceIndex,
    uint dimension) {
    const uvec4 block = raym0nadeCounterRandomBlockForDimension(
        seed, pixelIndex, sampleIndex, bounceIndex, dimension);
    return raym0nadeCounterRandomBlockWord(block, dimension);
}

float raym0nadeCounterRandomWordToOpen01(uint word) {
    return (float(word >> 9U) + 0.5) * (1.0 / 8388608.0);
}

float raym0nadeCounterRandomBlockOpen01(uvec4 block, uint dimension) {
    return raym0nadeCounterRandomWordToOpen01(
        raym0nadeCounterRandomBlockWord(block, dimension));
}

float raym0nadeCounterRandomOpen01(
    uint seed,
    uint pixelIndex,
    uint sampleIndex,
    uint bounceIndex,
    uint dimension) {
    return raym0nadeCounterRandomWordToOpen01(raym0nadeCounterRandomUint32(
        seed, pixelIndex, sampleIndex, bounceIndex, dimension));
}

#endif
