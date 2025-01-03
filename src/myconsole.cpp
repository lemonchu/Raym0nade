#include <sstream>
#include <iostream>
#include "myconsole.h"

MyConsole::MyConsole() = default;

void MyConsole::createModel(const std::string &model_id) {
    if (models.find(model_id) != models.end()) {
        std::cout << "Model (" << model_id << ") is already exists." << std::endl;
        return;
    }

    std::string model_folder, model_name, skyMap_name;
    std::cout << "Enter the model path (e.g., fbx/): ";
    std::cin >> model_folder;
    std::cout << "Enter the model name (e.g., model.fbx): ";
    std::cin >> model_name;
    std::cout << "Enter the sky map name (e.g., sky.hdr): ";
    std::cin >> skyMap_name;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    models.emplace(std::piecewise_construct,
                   std::forward_as_tuple(model_id),
                   std::forward_as_tuple(model_folder, model_name, skyMap_name));
    std::cout << "Model (" << model_id << ") created." << std::endl;
}

void MyConsole::createRenderArgs(const std::string &str) {
    if (renderArgs.find(str) != renderArgs.end()) {
        std::cout << "Args (" << str << ") is already exists." << std::endl;
        return;
    }

    vec3 direction, right, up;
    float D, R, U, accuracy, exposure, P_Direct;
    int width, height, spp, threads;
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

    std::cout << "exposure: ";
    std::cin >> exposure;

    std::cout << "width, height: ";
    std::cin >> width >> height;

    std::cout << "spp: ";
    std::cin >> spp;

    std::cout << "threads: ";
    std::cin >> threads;

    std::cout << "P_Direct: ";
    std::cin >> P_Direct;

    std::cout << "savePath: ";
    std::cin >> savePath;

    renderArgs[str] = {position, direction, up, right, accuracy, exposure, P_Direct, width, height, spp, threads, savePath};
    std::cout << "RenderArgs (" << str << ") created." << std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

void MyConsole::deleteModel(const std::string &str) {
    if (models.find(str) == models.end()) {
        std::cout << "Model (" << str << ") does not exists." << std::endl;
        return;
    }
    models.erase(str);
    std::cout << "Model (" << str << ") deleted." << std::endl;
}

void MyConsole::deleteRenderArgs(const std::string &str) {
    if (renderArgs.find(str) == renderArgs.end()) {
        std::cout << "Args (" << str << ") does not exists." << std::endl;
        return;
    }
    renderArgs.erase(str);
    std::cout << "RenderArgs (" << str << ") deleted." << std::endl;
}

void MyConsole::viewModel(const std::string &str) {
    if (models.find(str) == models.end()) {
        std::cout << "Model (" << str << ") does not exists." << std::endl;
        return;
    }
    std::cout << "Model Path: " << models[str].model_path << std::endl;
    std::cout << "Faces: " << models[str].faces.size() << std::endl;
}

void MyConsole::viewRenderArgs(const std::string &str) {
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
    std::cout << "exposure: " << args.exposure << std::endl;
    std::cout << "width, height: " << args.width << " " << args.height << std::endl;
    std::cout << "spp: " << args.spp << std::endl;
    std::cout << "threads: " << args.threads << std::endl;
    std::cout << "P_Direct: " << args.P_Direct << std::endl;
    std::cout << "savePath: " << args.savePath << std::endl;
}

void MyConsole::render(const std::string &model_str, const std::string &args_str) {
    if (models.find(model_str) == models.end()) {
        std::cout << "Model (" << model_str << ") does not exists." << std::endl;
        return;
    }
    if (renderArgs.find(args_str) == renderArgs.end()) {
        std::cout << "Args (" << args_str << ") does not exists." << std::endl;
        return;
    }
    render_multiThread(models[model_str], renderArgs[args_str]);
}

void parseCommand(MyConsole &console, const std::string &opt) {
    std::istringstream iss(opt);
    std::string command, type, name, model_str, args_str;
    iss >> command;

    if (command == "create") {
        iss >> type >> name;
        if (type == "model") {
            console.createModel(name);
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