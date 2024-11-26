#ifndef MYCONSOLE_H
#define MYCONSOLE_H

#include <map>
#include <string>
#include "model.h"
#include "render.h"

class MyConsole {
private:
    std::map<std::string, Model> models;
    std::map<std::string, RenderArgs> renderArgs;

public:
    MyConsole();

    void createModel(const std::string &model_id);

    void createRenderArgs(const std::string &str);

    void deleteModel(const std::string &str);

    void deleteRenderArgs(const std::string &str);

    void viewModel(const std::string &str);

    void viewRenderArgs(const std::string &str);

    void render(const std::string &model_str, const std::string &args_str);
};

void parseCommand(MyConsole &console, const std::string &opt);

#endif // MYCONSOLE_H