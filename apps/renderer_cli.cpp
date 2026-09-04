#include <iostream>

#include "raym0nade/console.hpp"

int main() {
    raym0nade::ConsoleApplication application{std::cin, std::cout, std::cerr};
    return application.run();
}
