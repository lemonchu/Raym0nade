#include "raym0nade/console.hpp"

#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace raym0nade {
namespace {

void requireIdentifier(const std::string& identifier, const char* kind) {
    if (identifier.empty()) {
        throw std::invalid_argument(std::string("Missing ") + kind + " identifier.");
    }
}

void requireInput(const std::istream& input, const char* field) {
    if (!input) {
        throw std::runtime_error(std::string("Failed to read ") + field + '.');
    }
}

}  // namespace

ConsoleApplication::ConsoleApplication(
    std::istream& input, std::ostream& output, std::ostream& error)
    : input_(input), output_(output), error_(error) {}

int ConsoleApplication::run() {
    std::string commandLine;
    while (true) {
        output_ << "> " << std::flush;
        if (!std::getline(input_, commandLine)) {
            return input_.eof() ? 0 : 1;
        }
        if (commandLine.empty()) {
            continue;
        }
        if (!execute(commandLine)) {
            return 0;
        }
    }
}

bool ConsoleApplication::execute(const std::string& commandLine) {
    try {
        std::istringstream command{commandLine};
        std::string verb;
        std::string kind;
        std::string firstId;
        std::string secondId;
        command >> verb;

        if (verb == "exit") {
            return false;
        }
        if (verb == "create") {
            command >> kind >> firstId;
            if (kind == "model") {
                createModel(firstId);
            } else if (kind == "settings" || kind == "args") {
                createRenderSettings(firstId);
            } else {
                throw std::invalid_argument("Usage: create <model|settings> <id>");
            }
            return true;
        }
        if (verb == "delete") {
            command >> kind >> firstId;
            if (kind == "model") {
                deleteModel(firstId);
            } else if (kind == "settings" || kind == "args") {
                deleteRenderSettings(firstId);
            } else {
                throw std::invalid_argument("Usage: delete <model|settings> <id>");
            }
            return true;
        }
        if (verb == "show" || verb == "view") {
            command >> kind >> firstId;
            if (kind == "model") {
                showModel(firstId);
            } else if (kind == "settings" || kind == "args") {
                showRenderSettings(firstId);
            } else {
                throw std::invalid_argument("Usage: show <model|settings> <id>");
            }
            return true;
        }
        if (verb == "render") {
            command >> firstId >> secondId;
            render(firstId, secondId);
            return true;
        }

        throw std::invalid_argument("Unknown command: " + verb);
    } catch (const std::exception& exception) {
        error_ << "Error: " << exception.what() << '\n';
        return true;
    }
}

void ConsoleApplication::createModel(const std::string& id) {
    requireIdentifier(id, "model");
    if (models_.find(id) != models_.end()) {
        throw std::invalid_argument("Model already exists: " + id);
    }

    std::filesystem::path directory;
    std::filesystem::path filename;
    std::filesystem::path skyFilename;
    output_ << "Model directory: " << std::flush;
    input_ >> directory;
    requireInput(input_, "model directory");
    output_ << "Model filename: " << std::flush;
    input_ >> filename;
    requireInput(input_, "model filename");
    output_ << "Sky filename (or null): " << std::flush;
    input_ >> skyFilename;
    requireInput(input_, "sky filename");
    input_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    auto model = std::make_unique<Model>(directory, filename, skyFilename);
    models_.emplace(id, std::move(model));
    output_ << "Model created: " << id << '\n';
}

void ConsoleApplication::createRenderSettings(const std::string& id) {
    requireIdentifier(id, "settings");
    if (renderSettings_.find(id) != renderSettings_.end()) {
        throw std::invalid_argument("Render settings already exist: " + id);
    }

    RenderSettings settings;
    const auto readVector = [this](vec3& value, const char* field) {
        input_ >> value.x >> value.y >> value.z;
        requireInput(input_, field);
    };

    output_ << "Direction (x y z): " << std::flush;
    readVector(settings.direction, "direction");
    output_ << "Right (x y z): " << std::flush;
    readVector(settings.right, "right vector");
    output_ << "Up (x y z): " << std::flush;
    readVector(settings.up, "up vector");

    output_ << "Position coefficients (direction right up): " << std::flush;
    float alongDirection = 0.0F;
    float alongRight = 0.0F;
    float alongUp = 0.0F;
    input_ >> alongDirection >> alongRight >> alongUp;
    requireInput(input_, "position coefficients");
    settings.position = alongDirection * settings.direction + alongRight * settings.right +
                        alongUp * settings.up;

    output_ << "Pixel scale, focus distance, circle of confusion, exposure: " << std::flush;
    input_ >> settings.pixelScale >> settings.focusDistance >> settings.circleOfConfusion >>
        settings.exposure;
    requireInput(input_, "camera scalars");
    output_ << "Width and height: " << std::flush;
    input_ >> settings.width >> settings.height;
    requireInput(input_, "image dimensions");
    output_ << "Samples per pixel, threads, direct-light probability: " << std::flush;
    input_ >> settings.samplesPerPixel >> settings.threadCount >> settings.directLightProbability;
    requireInput(input_, "sampling settings");
    output_ << "Output prefix: " << std::flush;
    input_ >> settings.outputPrefix;
    requireInput(input_, "output prefix");
    input_.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    settings.validate();
    renderSettings_.emplace(id, std::move(settings));
    output_ << "Render settings created: " << id << '\n';
}

void ConsoleApplication::deleteModel(const std::string& id) {
    requireIdentifier(id, "model");
    if (models_.erase(id) == 0U) {
        throw std::invalid_argument("Unknown model: " + id);
    }
    output_ << "Model deleted: " << id << '\n';
}

void ConsoleApplication::deleteRenderSettings(const std::string& id) {
    requireIdentifier(id, "settings");
    if (renderSettings_.erase(id) == 0U) {
        throw std::invalid_argument("Unknown render settings: " + id);
    }
    output_ << "Render settings deleted: " << id << '\n';
}

void ConsoleApplication::showModel(const std::string& id) const {
    requireIdentifier(id, "model");
    const auto model = models_.find(id);
    if (model == models_.end()) {
        throw std::invalid_argument("Unknown model: " + id);
    }
    output_ << "Model path: " << model->second->modelPath().string() << '\n';
    output_ << "Faces: " << model->second->faceCount() << '\n';
}

void ConsoleApplication::showRenderSettings(const std::string& id) const {
    requireIdentifier(id, "settings");
    const auto found = renderSettings_.find(id);
    if (found == renderSettings_.end()) {
        throw std::invalid_argument("Unknown render settings: " + id);
    }
    const RenderSettings& settings = found->second;
    output_ << "Direction: " << settings.direction.x << ' ' << settings.direction.y << ' '
            << settings.direction.z << '\n';
    output_ << "Right: " << settings.right.x << ' ' << settings.right.y << ' '
            << settings.right.z << '\n';
    output_ << "Up: " << settings.up.x << ' ' << settings.up.y << ' ' << settings.up.z << '\n';
    output_ << "Position: " << settings.position.x << ' ' << settings.position.y << ' '
            << settings.position.z << '\n';
    output_ << "Pixel scale: " << settings.pixelScale << '\n';
    output_ << "Exposure: " << settings.exposure << '\n';
    output_ << "Resolution: " << settings.width << 'x' << settings.height << '\n';
    output_ << "Samples per pixel: " << settings.samplesPerPixel << '\n';
    output_ << "Threads: " << settings.resolvedThreadCount() << '\n';
    output_ << "Direct-light probability: " << settings.directLightProbability << '\n';
    output_ << "Output prefix: " << settings.outputPrefix.string() << '\n';
}

void ConsoleApplication::render(const std::string& modelId, const std::string& settingsId) {
    requireIdentifier(modelId, "model");
    requireIdentifier(settingsId, "settings");
    const auto model = models_.find(modelId);
    const auto settings = renderSettings_.find(settingsId);
    if (model == models_.end()) {
        throw std::invalid_argument("Unknown model: " + modelId);
    }
    if (settings == renderSettings_.end()) {
        throw std::invalid_argument("Unknown render settings: " + settingsId);
    }
    (void)renderToFiles(*model->second, settings->second);
}

}  // namespace raym0nade
