cmake_minimum_required(VERSION 3.25)

project(
    Raym0nade
    VERSION 1.0.0
    DESCRIPTION "A portable offline path tracer"
    LANGUAGES CXX
)

# Conda and the other prebuilt Windows dependencies use the release dynamic
# runtime. Keep Debug optimization and symbols while using the same redistributable
# CRT ABI, so local and CI test binaries do not require the debug-only MSVC DLLs.
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
endif()

include(CTest)
include("${PROJECT_SOURCE_DIR}/cmake/CompilerWarnings.cmake")

option(RAYM0NADE_BUILD_RENDERER "Build the interactive renderer" ON)
option(RAYM0NADE_BUILD_FXAA "Build the standalone FXAA tool" ON)

find_package(assimp 5 CONFIG REQUIRED)
find_package(PNG REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Python is embedded temporarily for the HDR and DDS decoder bridge.
set(Python3_FIND_VIRTUALENV FIRST)
find_package(Python3 3.10 COMPONENTS Interpreter Development.Embed REQUIRED)

# Most Windows Python distributions, including Conda, do not ship a debug
# interpreter library. Python 3.12 and earlier select pythonXY_d.lib whenever
# _DEBUG is visible while Python.h is parsed. Force-include a small wrapper for
# the bridge translation unit so that Python's release ABI is selected without
# changing the debug ABI used by the rest of that source file.
if(MSVC AND NOT Python3_LIBRARY_DEBUG)
    set(RAYM0NADE_PYTHON_DEBUG_COMPATIBILITY_HEADER
        "${PROJECT_BINARY_DIR}/generated/raym0nade/python_debug_compatibility.hpp"
    )
    get_filename_component(
        RAYM0NADE_PYTHON_DEBUG_COMPATIBILITY_DIRECTORY
        "${RAYM0NADE_PYTHON_DEBUG_COMPATIBILITY_HEADER}"
        DIRECTORY
    )
    file(MAKE_DIRECTORY "${RAYM0NADE_PYTHON_DEBUG_COMPATIBILITY_DIRECTORY}")
    configure_file(
        "${PROJECT_SOURCE_DIR}/cmake/PythonDebugCompatibility.hpp.in"
        "${RAYM0NADE_PYTHON_DEBUG_COMPATIBILITY_HEADER}"
        COPYONLY
    )
    set_property(SOURCE src/python_image_loader.cpp APPEND PROPERTY
        COMPILE_OPTIONS
            "$<$<CONFIG:Debug>:/FI${RAYM0NADE_PYTHON_DEBUG_COMPATIBILITY_HEADER}>"
    )
endif()

set(RAYM0NADE_GLM_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/lib/glm")
if(NOT EXISTS "${RAYM0NADE_GLM_INCLUDE_DIR}/glm/glm.hpp")
    message(FATAL_ERROR
        "The GLM submodule is missing. Run: "
        "git submodule update --init --recursive"
    )
endif()

add_library(raym0nade_glm INTERFACE)
add_library(Raym0nade::glm ALIAS raym0nade_glm)
target_include_directories(raym0nade_glm SYSTEM INTERFACE
    "${RAYM0NADE_GLM_INCLUDE_DIR}"
)

add_library(raym0nade_core STATIC)
add_library(Raym0nade::core ALIAS raym0nade_core)

target_sources(raym0nade_core PRIVATE
    src/bvh.cpp
    src/component.cpp
    src/geometry.cpp
    src/image.cpp
    src/material.cpp
    src/model.cpp
    src/python_image_loader.cpp
    src/render.cpp
    src/sampling.cpp
)

target_compile_features(raym0nade_core PUBLIC cxx_std_17)
set_target_properties(raym0nade_core PROPERTIES
    CXX_EXTENSIONS OFF
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(raym0nade_core
    PUBLIC
        "${PROJECT_SOURCE_DIR}/include"
    PRIVATE
        "${PROJECT_SOURCE_DIR}/src"
)

file(TO_CMAKE_PATH "${PROJECT_SOURCE_DIR}/scripts" RAYM0NADE_SCRIPT_DIR)
target_compile_definitions(raym0nade_core PRIVATE
    "RAYM0NADE_SCRIPT_DIR=\"${RAYM0NADE_SCRIPT_DIR}\""
)

target_link_libraries(raym0nade_core
    PUBLIC
        Raym0nade::glm
    PRIVATE
        assimp::assimp
        PNG::PNG
        Python3::Python
        Threads::Threads
)

raym0nade_enable_warnings(raym0nade_core)

function(raym0nade_set_runtime_directory target)
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin"
    )

    foreach(configuration Debug Release RelWithDebInfo MinSizeRel)
        string(TOUPPER "${configuration}" configuration_upper)
        set_target_properties(${target} PROPERTIES
            "RUNTIME_OUTPUT_DIRECTORY_${configuration_upper}" "${PROJECT_BINARY_DIR}/bin"
        )
    endforeach()
endfunction()

if(RAYM0NADE_BUILD_RENDERER)
    add_executable(raym0nade_renderer
        apps/renderer_cli.cpp
        src/console.cpp
    )
    target_link_libraries(raym0nade_renderer PRIVATE Raym0nade::core)
    set_target_properties(raym0nade_renderer PROPERTIES OUTPUT_NAME raym0nade)
    raym0nade_enable_warnings(raym0nade_renderer)
    raym0nade_set_runtime_directory(raym0nade_renderer)
endif()

if(RAYM0NADE_BUILD_FXAA)
    add_executable(raym0nade_fxaa apps/fxaa_cli.cpp)
    target_link_libraries(raym0nade_fxaa PRIVATE Raym0nade::core)
    raym0nade_enable_warnings(raym0nade_fxaa)
    raym0nade_set_runtime_directory(raym0nade_fxaa)
endif()

if(BUILD_TESTING)
    add_executable(raym0nade_tests tests/core_tests.cpp)
    target_link_libraries(raym0nade_tests PRIVATE Raym0nade::core)
    raym0nade_enable_warnings(raym0nade_tests)
    raym0nade_set_runtime_directory(raym0nade_tests)

    add_executable(raym0nade_render_tests tests/render_smoke_tests.cpp)
    target_link_libraries(raym0nade_render_tests PRIVATE Raym0nade::core)
    target_compile_definitions(raym0nade_render_tests PRIVATE
        "RAYM0NADE_TEST_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\""
        "RAYM0NADE_TEST_BINARY_DIR=\"${PROJECT_BINARY_DIR}\""
    )
    raym0nade_enable_warnings(raym0nade_render_tests)
    raym0nade_set_runtime_directory(raym0nade_render_tests)

    add_executable(raym0nade_console_tests
        tests/console_tests.cpp
        src/console.cpp
    )
    target_link_libraries(raym0nade_console_tests PRIVATE Raym0nade::core)
    raym0nade_enable_warnings(raym0nade_console_tests)
    raym0nade_set_runtime_directory(raym0nade_console_tests)

    add_test(NAME core COMMAND raym0nade_tests)
    add_test(NAME render-smoke COMMAND raym0nade_render_tests)
    add_test(NAME console COMMAND raym0nade_console_tests)

    if(WIN32)
        get_filename_component(
            RAYM0NADE_PYTHON_RUNTIME_DIR
            "${Python3_EXECUTABLE}"
            DIRECTORY
        )
        set_tests_properties(core render-smoke console PROPERTIES
            ENVIRONMENT_MODIFICATION
                "PATH=path_list_prepend:${RAYM0NADE_PYTHON_RUNTIME_DIR};PATH=path_list_prepend:${RAYM0NADE_PYTHON_RUNTIME_DIR}/Library/bin"
        )
    endif()
endif()
