cmake_minimum_required(VERSION 3.10)
project(Raym0nade VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(SDL2 REQUIRED)
find_package(assimp CONFIG REQUIRED)
find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(libjpeg-turbo CONFIG REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

find_library(FREEIMAGE_LIBRARY freeimage REQUIRED)
find_path(FREEIMAGE_INCLUDE_DIR FreeImage.h REQUIRED)
include_directories(${FREEIMAGE_INCLUDE_DIR})

include_directories(
        include
        "${CMAKE_CURRENT_SOURCE_DIR}/glm"
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
        SDL2::SDL2
        assimp::assimp
        PNG::PNG
        ZLIB::ZLIB
        libjpeg-turbo::turbojpeg
        ${FREEIMAGE_LIBRARY}
        Python3::Python
)