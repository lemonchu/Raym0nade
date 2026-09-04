#include <iostream>
#include <sstream>
#include <string_view>

#include "raym0nade/console.hpp"

namespace {

int failureCount = 0;

void expect(bool condition, std::string_view message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        ++failureCount;
    }
}

void testWhitespaceAroundExit() {
    std::istringstream input{"  exit  \n"};
    std::ostringstream output;
    std::ostringstream error;
    raym0nade::ConsoleApplication application{input, output, error};

    expect(
        application.run() == 0,
        "A whitespace-padded exit command must terminate successfully.");
    expect(error.str().empty(), "A whitespace-padded exit command must not report an error.");
}

}  // namespace

int main() {
    testWhitespaceAroundExit();

    if (failureCount != 0) {
        std::cerr << failureCount << " console test assertion(s) failed.\n";
        return 1;
    }

    std::cout << "All console tests passed.\n";
    return 0;
}
