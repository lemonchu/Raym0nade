cmake_minimum_required(VERSION 3.10)
project(Raym0nade VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(SDL2 REQUIRED)
find_package(assimp CONFIG REQUIRED)
find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(libjpeg-turbo CONFIG REQUIRED)

include_directories(
        include
        "${CMAKE_CURRENT_SOURCE_DIR}/glm"
)

file(GLOB_RECURSE PROJECT_SOURCES "src/*.cpp")

add_executable(Demo src/main.cpp ${PROJECT_SOURCES})

target_link_libraries(Demo PUBLIC
        SDL2::SDL2
        assimp::assimp
        PNG::PNG
        ZLIB::ZLIB
        libjpeg-turbo::turbojpeg
)