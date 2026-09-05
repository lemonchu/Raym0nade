cmake_minimum_required(VERSION 3.25)

project(
    Raym0nade
    VERSION 1.0.0
    DESCRIPTION "A portable offline path tracer"
    LANGUAGES CXX
)

# Some Simplified Chinese MSVC installations emit UTF-8 /showIncludes text that CMake
# transcodes a second time. Ninja then records zero header dependencies and can reuse ABI-stale
# objects. Repair that one detected byte pattern while leaving every other locale untouched.
if(MSVC AND CMAKE_GENERATOR MATCHES "Ninja")
    string(HEX "${CMAKE_CL_SHOWINCLUDES_PREFIX}" RAYM0NADE_SHOWINCLUDES_PREFIX_HEX)
    if(RAYM0NADE_SHOWINCLUDES_PREFIX_HEX STREQUAL
       "e5a889e384a6e589b03a20e98d96e591ade68388e98f82e59ba6e6aca23a2020")
        string(ASCII
            230 179 168 230 132 143 58 32
            229 140 133 229 144 171 230 150 135 228 187 182 58 32 32
            CMAKE_CL_SHOWINCLUDES_PREFIX
        )
    endif()
endif()

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
option(RAYM0NADE_ENABLE_VULKAN "Build the experimental Vulkan GPU backend" OFF)

find_package(assimp 5 CONFIG REQUIRED)
find_package(PNG REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if(RAYM0NADE_ENABLE_VULKAN)
    if(WIN32 AND DEFINED ENV{CONDA_PREFIX} AND NOT DEFINED Vulkan_ROOT)
        set(Vulkan_ROOT "$ENV{CONDA_PREFIX}/Library")
    endif()
    find_package(Vulkan 1.2.162 REQUIRED COMPONENTS glslc)
    find_program(
        RAYM0NADE_SPIRV_VAL_EXECUTABLE
        NAMES spirv-val
        HINTS
            "${Vulkan_ROOT}/bin"
            "$ENV{CONDA_PREFIX}/bin"
            "$ENV{CONDA_PREFIX}/Library/bin"
        REQUIRED
    )
endif()

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
    src/render_contract.cpp
    src/sampling.cpp
    src/scene_data.cpp
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

if(RAYM0NADE_ENABLE_VULKAN)
    set(RAYM0NADE_SHADER_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/generated/shaders")

    function(raym0nade_compile_vulkan_shader output_variable source_filename description)
        set(source_path "${PROJECT_SOURCE_DIR}/shaders/${source_filename}")
        set(output_path "${RAYM0NADE_SHADER_OUTPUT_DIRECTORY}/${source_filename}.spv")
        add_custom_command(
            OUTPUT "${output_path}"
            COMMAND "${CMAKE_COMMAND}" -E make_directory "${RAYM0NADE_SHADER_OUTPUT_DIRECTORY}"
            COMMAND Vulkan::glslc
                --target-env=vulkan1.2
                "${source_path}"
                -o "${output_path}"
            COMMAND "${RAYM0NADE_SPIRV_VAL_EXECUTABLE}"
                --target-env vulkan1.2
                "${output_path}"
            DEPENDS
                "${source_path}"
                Vulkan::glslc
            COMMENT "${description}"
            VERBATIM
        )
        set(${output_variable} "${output_path}" PARENT_SCOPE)
    endfunction()

    raym0nade_compile_vulkan_shader(
        RAYM0NADE_RAY_QUERY_SHADER_BINARY
        ray_query_triangle.comp
        "Compiling the Vulkan Ray Query self-test shader"
    )
    raym0nade_compile_vulkan_shader(
        RAYM0NADE_RAY_QUERY_SCENE_SHADER_BINARY
        ray_query_scene.comp
        "Compiling the reusable Vulkan Ray Query scene shader"
    )
    raym0nade_compile_vulkan_shader(
        RAYM0NADE_PRIMARY_AOV_SHADER_BINARY
        primary_aov.comp
        "Compiling the Vulkan primary-AOV shader"
    )
    add_custom_target(raym0nade_gpu_shaders
        DEPENDS
            "${RAYM0NADE_RAY_QUERY_SHADER_BINARY}"
            "${RAYM0NADE_RAY_QUERY_SCENE_SHADER_BINARY}"
            "${RAYM0NADE_PRIMARY_AOV_SHADER_BINARY}"
    )

    add_library(raym0nade_vulkan STATIC
        src/gpu/vulkan_capabilities.cpp
        src/gpu/vulkan_primary_renderer.cpp
        src/gpu/vulkan_ray_query.cpp
        src/gpu/vulkan_runtime.cpp
        src/gpu/vulkan_self_test.cpp
    )
    add_library(Raym0nade::vulkan ALIAS raym0nade_vulkan)
    target_compile_features(raym0nade_vulkan PUBLIC cxx_std_17)
    set_target_properties(raym0nade_vulkan PROPERTIES
        CXX_EXTENSIONS OFF
        POSITION_INDEPENDENT_CODE ON
    )
    target_include_directories(raym0nade_vulkan PUBLIC
        "${PROJECT_SOURCE_DIR}/include"
    )
    target_link_libraries(raym0nade_vulkan
        PUBLIC
            Raym0nade::core
        PRIVATE
            Vulkan::Vulkan
    )
    raym0nade_enable_warnings(raym0nade_vulkan)

    add_executable(raym0nade_gpu_probe apps/gpu_probe.cpp)
    target_link_libraries(raym0nade_gpu_probe PRIVATE Raym0nade::vulkan)
    add_dependencies(raym0nade_gpu_probe raym0nade_gpu_shaders)
    file(TO_CMAKE_PATH
        "${RAYM0NADE_RAY_QUERY_SHADER_BINARY}"
        RAYM0NADE_RAY_QUERY_SHADER_BINARY_CMAKE_PATH
    )
    target_compile_definitions(raym0nade_gpu_probe PRIVATE
        "RAYM0NADE_RAY_QUERY_SELF_TEST_SHADER=\"${RAYM0NADE_RAY_QUERY_SHADER_BINARY_CMAKE_PATH}\""
    )
    raym0nade_enable_warnings(raym0nade_gpu_probe)
    raym0nade_set_runtime_directory(raym0nade_gpu_probe)
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

    add_executable(raym0nade_primary_render_tests tests/primary_render_tests.cpp)
    target_link_libraries(raym0nade_primary_render_tests PRIVATE Raym0nade::core)
    target_compile_definitions(raym0nade_primary_render_tests PRIVATE
        "RAYM0NADE_TEST_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\""
    )
    raym0nade_enable_warnings(raym0nade_primary_render_tests)
    raym0nade_set_runtime_directory(raym0nade_primary_render_tests)

    add_executable(raym0nade_console_tests
        tests/console_tests.cpp
        src/console.cpp
    )
    target_link_libraries(raym0nade_console_tests PRIVATE Raym0nade::core)
    raym0nade_enable_warnings(raym0nade_console_tests)
    raym0nade_set_runtime_directory(raym0nade_console_tests)

    if(RAYM0NADE_ENABLE_VULKAN)
        add_executable(raym0nade_gpu_scene_tests tests/gpu_scene_tests.cpp)
        target_link_libraries(raym0nade_gpu_scene_tests PRIVATE
            Raym0nade::core
            Raym0nade::vulkan
        )
        add_dependencies(raym0nade_gpu_scene_tests raym0nade_gpu_shaders)
        file(TO_CMAKE_PATH
            "${RAYM0NADE_RAY_QUERY_SCENE_SHADER_BINARY}"
            RAYM0NADE_RAY_QUERY_SCENE_SHADER_BINARY_CMAKE_PATH
        )
        target_compile_definitions(raym0nade_gpu_scene_tests PRIVATE
            "RAYM0NADE_TEST_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\""
            "RAYM0NADE_RAY_QUERY_SCENE_SHADER=\"${RAYM0NADE_RAY_QUERY_SCENE_SHADER_BINARY_CMAKE_PATH}\""
        )
        raym0nade_enable_warnings(raym0nade_gpu_scene_tests)
        raym0nade_set_runtime_directory(raym0nade_gpu_scene_tests)

        add_executable(
            raym0nade_gpu_primary_render_tests
            tests/gpu_primary_render_tests.cpp
        )
        target_link_libraries(raym0nade_gpu_primary_render_tests PRIVATE
            Raym0nade::core
            Raym0nade::vulkan
        )
        add_dependencies(raym0nade_gpu_primary_render_tests raym0nade_gpu_shaders)
        file(TO_CMAKE_PATH
            "${RAYM0NADE_PRIMARY_AOV_SHADER_BINARY}"
            RAYM0NADE_PRIMARY_AOV_SHADER_BINARY_CMAKE_PATH
        )
        target_compile_definitions(raym0nade_gpu_primary_render_tests PRIVATE
            "RAYM0NADE_TEST_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\""
            "RAYM0NADE_PRIMARY_AOV_SHADER=\"${RAYM0NADE_PRIMARY_AOV_SHADER_BINARY_CMAKE_PATH}\""
        )
        raym0nade_enable_warnings(raym0nade_gpu_primary_render_tests)
        raym0nade_set_runtime_directory(raym0nade_gpu_primary_render_tests)
    endif()

    add_test(NAME core COMMAND raym0nade_tests)
    add_test(NAME render-smoke COMMAND raym0nade_render_tests)
    add_test(NAME primary-render COMMAND raym0nade_primary_render_tests)
    add_test(NAME console COMMAND raym0nade_console_tests)
    if(RAYM0NADE_ENABLE_VULKAN)
        add_test(NAME gpu-ray-query COMMAND raym0nade_gpu_probe --self-test)
        add_test(NAME gpu-packed-scene COMMAND raym0nade_gpu_scene_tests)
        add_test(NAME gpu-primary-render COMMAND raym0nade_gpu_primary_render_tests)
        set_tests_properties(gpu-ray-query gpu-packed-scene gpu-primary-render PROPERTIES
            SKIP_RETURN_CODE 77
        )
    endif()

    if(WIN32)
        get_filename_component(
            RAYM0NADE_PYTHON_RUNTIME_DIR
            "${Python3_EXECUTABLE}"
            DIRECTORY
        )
        set(RAYM0NADE_WINDOWS_RUNTIME_TESTS core render-smoke primary-render console)
        if(RAYM0NADE_ENABLE_VULKAN)
            list(APPEND RAYM0NADE_WINDOWS_RUNTIME_TESTS
                gpu-ray-query
                gpu-packed-scene
                gpu-primary-render
            )
        endif()
        set_tests_properties(${RAYM0NADE_WINDOWS_RUNTIME_TESTS} PROPERTIES
            ENVIRONMENT_MODIFICATION
                "PATH=path_list_prepend:${RAYM0NADE_PYTHON_RUNTIME_DIR};PATH=path_list_prepend:${RAYM0NADE_PYTHON_RUNTIME_DIR}/Library/bin"
        )
        if(RAYM0NADE_ENABLE_VULKAN)
            set_property(TEST gpu-ray-query gpu-packed-scene gpu-primary-render APPEND PROPERTY
                ENVIRONMENT_MODIFICATION
                    "VK_LAYER_PATH=set:${RAYM0NADE_PYTHON_RUNTIME_DIR}/Library/bin"
            )
        endif()
    endif()
endif()
