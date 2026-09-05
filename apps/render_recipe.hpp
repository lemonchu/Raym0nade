#pragma once

#include <filesystem>
#include <fstream>
#include <istream>
#include <stdexcept>
#include <string>

#include "raym0nade/render.hpp"

namespace raym0nade::cli {

struct SingleRenderRecipe {
    std::filesystem::path modelDirectory;
    std::filesystem::path modelFilename;
    std::filesystem::path skyFilename;
    RenderSettings settings;
};

[[nodiscard]] inline std::filesystem::path resolveFromWorkingDirectory(
    const std::filesystem::path& path,
    const std::filesystem::path& workingDirectory) {
    if (path.is_absolute()) {
        return path.lexically_normal();
    }
    return (workingDirectory / path).lexically_normal();
}

namespace detail {

inline void requireRecipeToken(
    std::istream& input,
    const std::string& expected,
    const char* context) {
    std::string actual;
    if (!(input >> actual) || actual != expected) {
        throw std::runtime_error(
            std::string{"Invalid recipe while reading "} + context +
            ": expected '" + expected + "'.");
    }
}

template <typename Value>
void requireRecipeValue(
    std::istream& input,
    Value& value,
    const char* description) {
    if (!(input >> value)) {
        throw std::runtime_error(
            std::string{"Invalid recipe value for "} + description + '.');
    }
}

inline void readRecipeVector(
    std::istream& input,
    vec3& value,
    const char* description) {
    requireRecipeValue(input, value.x, description);
    requireRecipeValue(input, value.y, description);
    requireRecipeValue(input, value.z, description);
}

}  // namespace detail

[[nodiscard]] inline SingleRenderRecipe loadSingleRenderRecipe(
    const std::filesystem::path& filename,
    const std::filesystem::path& workingDirectory) {
    std::ifstream input{filename};
    if (!input) {
        throw std::runtime_error(
            "Could not open render recipe: " + filename.string());
    }

    SingleRenderRecipe result;
    std::string modelId;
    std::string settingsId;
    detail::requireRecipeToken(input, "create", "model command");
    detail::requireRecipeToken(input, "model", "model command");
    detail::requireRecipeValue(input, modelId, "model identifier");
    detail::requireRecipeValue(input, result.modelDirectory, "model directory");
    detail::requireRecipeValue(input, result.modelFilename, "model filename");
    detail::requireRecipeValue(input, result.skyFilename, "sky filename");

    detail::requireRecipeToken(input, "create", "settings command");
    std::string settingsCommand;
    detail::requireRecipeValue(input, settingsCommand, "settings command type");
    if (settingsCommand != "settings" && settingsCommand != "args") {
        throw std::runtime_error(
            "Invalid recipe while reading settings command: expected 'settings' or 'args'.");
    }
    detail::requireRecipeValue(input, settingsId, "settings identifier");
    detail::readRecipeVector(input, result.settings.direction, "camera direction");
    detail::readRecipeVector(input, result.settings.right, "camera right vector");
    detail::readRecipeVector(input, result.settings.up, "camera up vector");

    float alongDirection = 0.0F;
    float alongRight = 0.0F;
    float alongUp = 0.0F;
    detail::requireRecipeValue(
        input,
        alongDirection,
        "camera position direction coefficient");
    detail::requireRecipeValue(
        input,
        alongRight,
        "camera position right coefficient");
    detail::requireRecipeValue(
        input,
        alongUp,
        "camera position up coefficient");
    result.settings.position =
        alongDirection * result.settings.direction +
        alongRight * result.settings.right +
        alongUp * result.settings.up;

    detail::requireRecipeValue(input, result.settings.pixelScale, "pixel scale");
    detail::requireRecipeValue(input, result.settings.focusDistance, "focus distance");
    detail::requireRecipeValue(
        input,
        result.settings.circleOfConfusion,
        "circle of confusion");
    detail::requireRecipeValue(input, result.settings.exposure, "exposure");
    detail::requireRecipeValue(input, result.settings.width, "width");
    detail::requireRecipeValue(input, result.settings.height, "height");
    detail::requireRecipeValue(
        input,
        result.settings.samplesPerPixel,
        "samples per pixel");
    detail::requireRecipeValue(
        input,
        result.settings.threadCount,
        "CPU thread count");
    detail::requireRecipeValue(
        input,
        result.settings.directLightProbability,
        "direct-light probability");
    detail::requireRecipeValue(
        input,
        result.settings.outputPrefix,
        "output prefix");

    detail::requireRecipeToken(input, "render", "render command");
    std::string renderedModelId;
    std::string renderedSettingsId;
    detail::requireRecipeValue(
        input,
        renderedModelId,
        "render model identifier");
    detail::requireRecipeValue(
        input,
        renderedSettingsId,
        "render settings identifier");
    if (renderedModelId != modelId || renderedSettingsId != settingsId) {
        throw std::runtime_error(
            "The recipe render command does not reference its created objects.");
    }

    std::string trailing;
    if (input >> trailing) {
        if (trailing != "exit") {
            throw std::runtime_error(
                "Unexpected trailing recipe token: " + trailing);
        }
        if (input >> trailing) {
            throw std::runtime_error(
                "Unexpected content after the recipe exit command.");
        }
    }

    result.modelDirectory =
        resolveFromWorkingDirectory(result.modelDirectory, workingDirectory);
    result.modelFilename =
        result.modelFilename.is_absolute()
            ? result.modelFilename.lexically_normal()
            : (result.modelDirectory / result.modelFilename).lexically_normal();
    if (!result.skyFilename.empty() && result.skyFilename != "null") {
        result.skyFilename =
            result.skyFilename.is_absolute()
                ? result.skyFilename.lexically_normal()
                : (result.modelFilename.parent_path() / result.skyFilename)
                      .lexically_normal();
    }
    result.settings.outputPrefix = resolveFromWorkingDirectory(
        result.settings.outputPrefix,
        workingDirectory);
    return result;
}

}  // namespace raym0nade::cli
