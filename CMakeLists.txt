cmake_minimum_required(VERSION 3.10)
project(Raym0nade VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(SDL2 REQUIRED)

message(STATUS "SDL2_INCLUDE_DIRS: ${SDL2_INCLUDE_DIRS}")
message(STATUS "SDL2_LIBRARIES: ${SDL2_LIBRARIES}")

include_directories(${SDL2_INCLUDE_DIRS})
include_directories(include)

find_package(assimp CONFIG REQUIRED)
find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

file(GLOB_RECURSE PROJECT_SOURCES "src/*.cpp" "src/geometry/*.cpp" "src/mathematics/*.cpp")

message(STATUS "PROJECT_SOURCES: ${PROJECT_SOURCES}")

add_executable(Demo src/main.cpp ${PROJECT_SOURCES})

target_link_libraries(Demo PUBLIC ${SDL2_LIBRARIES})
target_link_libraries(Demo PUBLIC assimp::assimp)
target_link_libraries(Demo PUBLIC PNG::PNG)
target_link_libraries(Demo PUBLIC ${ZLIB_LIBRARIES})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/glm")
