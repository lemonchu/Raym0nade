cmake_minimum_required(VERSION 3.10)
project(Raym0nade VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(assimp CONFIG REQUIRED)
find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

include_directories(
        include
        "${CMAKE_CURRENT_SOURCE_DIR}/lib/glm"
        ${Python3_INCLUDE_DIRS}
)

file(GLOB_RECURSE PROJECT_SOURCES "src/*.cpp")

add_executable(Demo main.cpp ${PROJECT_SOURCES})

# Set the working directory to the root directory
set_target_properties(Demo PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_CURRENT_SOURCE_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(Demo PUBLIC
        assimp::assimp
        PNG::PNG
        ZLIB::ZLIB
        Python3::Python
)