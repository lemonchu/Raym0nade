#ifndef MYCONSOLE_H
#define MYCONSOLE_H

#include <sstream>
#include <iostream>
#include <map>
#include <string>
#include "model.h"
#include "render.h"

class MyConsole {
private:
    std::map<std::string, Model> models;
    std::map<std::string, RenderArgs> renderArgs;

public:
    void createModel(const std::string &model_id, const std::string &model_path, const std::string &model_name) {
        if (models.find(model_id) != models.end()) {
            std::cout << "Model (" << model_id << ") is already exists." << std::endl;
            return;
        }
        models[model_id] = Model();
        if (models[model_id].load(model_path, model_name) != 0) {
            std::cout << "Failed to load model (" << model_id << ")." << std::endl;
            models.erase(model_id);
        } else {
            std::cout << "Model (" << model_id << ") created." << std::endl;
        }
    }

    void createRenderArgs(const std::string &str) {
        if (renderArgs.find(str) != renderArgs.end()) {
            std::cout << "Args (" << str << ") is already exists." << std::endl;
            return;
        }

        vec3 direction, right, up;
        float D, R, U, accuracy;
        unsigned int oversampling, width, height;
        std::string savePath;

        std::cout << "direction (x,y,z): ";
        std::cin >> direction.x >> direction.y >> direction.z;

        std::cout << "right (x,y,z): ";
        std::cin >> right.x >> right.y >> right.z;

        std::cout << "up (x,y,z): ";
        std::cin >> up.x >> up.y >> up.z;

        std::cout << "position (D,R,U): ";
        std::cin >> D >> R >> U;
        glm::vec3 position = D * direction + R * right + U * up;

        std::cout << "accuracy: ";
        std::cin >> accuracy;

        std::cout << "Oversampling: ";
        std::cin >> oversampling;

        std::cout << "width, height: ";
        std::cin >> width >> height;

        std::cout << "savePath: ";
        std::cin >> savePath;

        renderArgs[str] = {position, direction, up, right, accuracy, oversampling, width, height, savePath};
        std::cout << "RenderArgs (" << str << ") created." << std::endl;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    void deleteModel(const std::string &str) {
        if (models.find(str) == models.end()) {
            std::cout << "Model (" << str << ") does not exists." << std::endl;
            return;
        }
        models.erase(str);
        std::cout << "Model (" << str << ") deleted." << std::endl;
    }

    void deleteRenderArgs(const std::string &str) {
        if (renderArgs.find(str) == renderArgs.end()) {
            std::cout << "Args (" << str << ") does not exists." << std::endl;
            return;
        }
        renderArgs.erase(str);
        std::cout << "RenderArgs (" << str << ") deleted." << std::endl;
    }

    void viewModel(const std::string &str) {
        if (models.find(str) == models.end()) {
            std::cout << "Model (" << str << ") does not exists." << std::endl;
            return;
        }
        std::cout << "Model path: " << models[str].model_path << std::endl;
        std::cout << "Faces: " << models[str].triangles.size() << std::endl;
    }

    void viewRenderArgs(const std::string &str) {
        if (renderArgs.find(str) == renderArgs.end()) {
            std::cout << "Args (" << str << ") does not exists." << std::endl;
            return;
        }
        const RenderArgs &args = renderArgs[str];
        std::cout << "direction : " << args.direction.x << " " << args.direction.y << " " << args.direction.z
                  << std::endl;
        std::cout << "right : " << args.right.x << " " << args.right.y << " " << args.right.z << std::endl;
        std::cout << "up : " << args.up.x << " " << args.up.y << " " << args.up.z << std::endl;
        std::cout << "position : " << args.position.x << " " << args.position.y << " " << args.position.z << std::endl;
        std::cout << "accuracy: " << args.accuracy << std::endl;
        std::cout << "Oversampling: " << args.oversampling << std::endl;
        std::cout << "width, height: " << args.width << " " << args.height << std::endl;
    }

    void render(const std::string &model_str, const std::string &args_str) {
        if (models.find(model_str) == models.end()) {
            std::cout << "Model (" << model_str << ") does not exists." << std::endl;
            return;
        }
        if (renderArgs.find(args_str) == renderArgs.end()) {
            std::cout << "Args (" << args_str << ") does not exists." << std::endl;
            return;
        }
        Renderer renderer;
        renderer.render(models[model_str], renderArgs[args_str]);
        std::cout << "Rendering completed." << std::endl;
    }
};

void parseCommand(MyConsole &console, const std::string &opt) {
    std::istringstream iss(opt);
    std::string command, type, name, model_str, args_str;
    iss >> command;

    if (command == "create") {
        iss >> type >> name;
        if (type == "model") {
            std::string model_path, model_name;
            std::cout << "Enter model path: ";
            std::cin >> model_path;
            std::cout << "Enter model name: ";
            std::cin >> model_name;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            console.createModel(name, model_path, model_name);
        } else if (type == "args") {
            console.createRenderArgs(name);
        }
    } else if (command == "delete") {
        iss >> type >> name;
        if (type == "model") {
            console.deleteModel(name);
        } else if (type == "args") {
            console.deleteRenderArgs(name);
        }
    } else if (command == "view") {
        iss >> type >> name;
        if (type == "model") {
            console.viewModel(name);
        } else if (type == "args") {
            console.viewRenderArgs(name);
        }
    } else if (command == "render") {
        iss >> model_str >> args_str;
        console.render(model_str, args_str);
    } else {
        std::cout << "Unknown command." << std::endl;
    }
}

#endif // MYCONSOLE_H