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

    RenderArgs &arg = renderArgs[str];

    auto read = [](vec3 &v) {
        std::cin >> v.x >> v.y >> v.z;
    };

    std::cout << "direction (x,y,z): ";
    read(arg.direction);

    std::cout << "right (x,y,z): ";
    read(arg.right);

    std::cout << "up (x,y,z): ";
    read(arg.up);

    std::cout << "position (D,R,U): ";
    float D, R, U;
    std::cin >> D >> R >> U;
    arg.position = D * arg.direction + R * arg.right + U * arg.up;

    std::cout << "accuracy, focus, CoC, exposure: ";
    std::cin >> arg.accuracy >> arg.focus >> arg.CoC >> arg.exposure;

    std::cout << "width, height: ";
    std::cin >> arg.width >> arg.height;

    std::cout << "spp, threads, P_Direct: : ";
    std::cin >> arg.spp >> arg.threads >> arg.P_Direct;

    std::cout << "savePath: ";
    std::cin >> arg.savePath;

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