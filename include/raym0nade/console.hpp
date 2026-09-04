#pragma once

#include <iosfwd>
#include <map>
#include <memory>
#include <string>

#include "raym0nade/render.hpp"

namespace raym0nade {

class ConsoleApplication {
public:
    ConsoleApplication(std::istream& input, std::ostream& output, std::ostream& error);

    int run();
    bool execute(const std::string& commandLine);

private:
    void createModel(const std::string& id);
    void createRenderSettings(const std::string& id);
    void deleteModel(const std::string& id);
    void deleteRenderSettings(const std::string& id);
    void showModel(const std::string& id) const;
    void showRenderSettings(const std::string& id) const;
    void render(const std::string& modelId, const std::string& settingsId);

    std::istream& input_;
    std::ostream& output_;
    std::ostream& error_;
    std::map<std::string, std::unique_ptr<Model>> models_;
    std::map<std::string, RenderSettings> renderSettings_;
};

}  // namespace raym0nade
