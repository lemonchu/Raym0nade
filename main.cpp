#include <iostream>
#include "render.h"
#include "myconsole.h"
#include <FreeImage.h>

int main() {
#if defined(__OPTIMIZE__) && !defined(__OPTIMIZE_SIZE__)
    std::cout << "O2 optimization is enabled." << std::endl;
#else
    std::cout << "O2 optimization is not enabled." << std::endl;
#endif
    MyConsole console;
    std::string opt;
    FreeImage_Initialise();
    while (true) {
        std::cout << "> ";
        std::getline(std::cin, opt);
        if (opt == "exit") {
            break;
        }
        if (opt == "") {
            continue;
        }
        parseCommand(console, opt);
    }
    FreeImage_DeInitialise();
    return 0;
}
