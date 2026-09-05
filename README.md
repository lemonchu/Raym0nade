# Raym0nade

Raym0nade is an offline path tracer written in C++17. Its CPU renderer is the correctness
reference, and its optional AMD Vulkan compute renderer executes the multi-bounce path estimator
with hardware Ray Queries. The project imports scenes with Assimp and provides direct and indirect
light sampling, textured materials, HDR environment lighting, denoising, bloom, depth of field,
and FXAA.

The current refactor establishes a portable core library and separate command-line applications.
Windows, Linux, and macOS use the same CMake presets and the same Conda dependency definition.

## Repository layout

```text
apps/                 Command-line entry points
cmake/                Reusable CMake helpers
include/raym0nade/    Public C++ headers
scripts/              Temporary Python image-decoder modules
src/                  Core implementation
tests/                Lightweight regression tests
```

The build creates these targets:

- `Raym0nade::core`: the reusable static renderer library
- `raym0nade_renderer`: the interactive renderer, output as `raym0nade`
- `raym0nade_fxaa`: the standalone FXAA utility
- `raym0nade_tests`: the CTest regression executable when testing is enabled
- `raym0nade_scene_data_tests`: packed-scene ABI, validation, texture deduplication, and mip tests
- `raym0nade_counter_rng_tests`: CPU counter-RNG known-answer and scheduling-order tests
- `raym0nade_render_tests`: the tiny lit imported-scene render regression executable when testing
  is enabled
- `raym0nade_primary_render_tests`: the deterministic, no-file CPU primary-AOV regression
  executable when testing is enabled
- `raym0nade_console_tests`: the command-loop regression executable when testing is enabled
- `raym0nade_vulkan`: the optional AMD Vulkan backend library
- `raym0nade_gpu_probe`: the optional AMD Ray Query capability and hardware self-test executable
- `raym0nade_gpu_scene_tests`: the optional packed-scene CPU/GPU intersection comparison
- `raym0nade_gpu_primary_render_tests`: the optional CPU/GPU primary-AOV comparison
- `raym0nade_gpu_counter_rng_tests`: the optional CPU/GLSL Philox known-answer comparison on a
  Vulkan device
- `raym0nade_gpu_path_render_tests`: the optional GPU path integration, determinism, partitioning,
  lighting, material, Film, and CPU-oracle regression executable
- `raym0nade_gpu_render`: the optional general AMD Vulkan path-rendering command-line application
- `raym0nade_gpu_primary_benchmark`: the optional Bistro ShapeNormal comparison and GPU-only
  diagnostic utility

See `docs/architecture.md` for module boundaries, ownership rules, render data flow,
determinism guarantees, current limitations, and the staged GPU roadmap. The AMD GPU backend is
specified in `docs/gpu-backend.md`.

## Requirements

All platforms require Git, Miniforge, and a native C++17 compiler:

- Windows: Visual Studio 2022 Build Tools with the Desktop development with C++ workload
- Linux: GCC 12 or newer, or a recent Clang toolchain
- macOS: Xcode Command Line Tools with Apple Clang

Clone the GLM submodule together with the repository:

```sh
git clone --recurse-submodules <repository-url>
```

For an existing checkout, initialize it with:

```sh
git submodule update --init --recursive
```

## Create the development environment

Before invoking Conda, clear any process-local `HTTP_PROXY`, `HTTPS_PROXY`, and `ALL_PROXY`
variables that route through a local proxy such as `127.0.0.1:7897`, and set `NO_PROXY=*` for that
shell. Do not change persistent user or system proxy settings. Then create and activate the
environment directly:

```sh
conda env create --file environment.yml
conda activate raym0nade
```

Update an existing environment with:

```sh
conda env update --name raym0nade --file environment.yml --prune
```

Assimp, libpng, Python, ImageIO, Pillow, NumPy, and FreeImage are isolated inside this environment.
The native and Python image packages are temporary dependencies while DDS and HDR decoding is
migrated to a native C++ implementation.

## Configure and build

Activate the `raym0nade` environment first. On Windows, run these commands from an x64 Native
Tools Command Prompt for Visual Studio 2022 so Ninja can find MSVC. The same CMake commands work on
all supported platforms:

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

The workflow preset combines all three steps:

```sh
cmake --workflow --preset release
```

On Windows, Debug builds keep unoptimized code, assertions, and debug symbols while using the
release dynamic MSVC runtime (`/MD`). This matches the ABI of the prebuilt Conda dependencies and
does not require Visual Studio's debug-only runtime DLLs.

Build products stay in a dedicated build tree:

```text
build/release/bin/raym0nade
build/release/bin/raym0nade_fxaa
build/release/bin/raym0nade_tests
build/release/bin/raym0nade_render_tests
build/release/bin/raym0nade_primary_render_tests
build/release/bin/raym0nade_console_tests
```

Windows appends `.exe` to these file names. Keep the Conda environment active while running a
program so its dynamic dependencies can be located.

## Run the renderer

With the Conda environment active, run the executable directly. On Windows:

```bat
build\release\bin\raym0nade.exe
```

On Linux or macOS:

```sh
./build/release/bin/raym0nade
```

The interactive command sequence is:

```text
create model <model-name>
create settings <settings-name>
render <model-name> <settings-name>
exit
```

See `docs/render-examples.md` for validated current recipes, output-selection guidance, and clearly
separated archival legacy settings. Ready-to-pipe command input is stored under `examples/`. Large
models and generated images belong under `model/` and `output/`; both directories are intentionally
excluded from version control.

## Tests and continuous integration

CPU CTest presets run the core regression suite, the tiny lit cross-thread imported-scene render
smoke test, the deterministic no-file primary-AOV oracle, and the command-loop regression suite:

```sh
ctest --preset debug
ctest --preset release
```

The GitHub Actions workflow configures, builds, and tests Release on Windows, Ubuntu, and macOS. It
also builds and tests Windows Debug to cover the Conda-compatible `/MD` ABI configuration. No
dependency is downloaded by CMake itself; external libraries are resolved from the active Conda
environment, and GLM comes from the pinned Git submodule.

## Known import limitation

Assimp 5.4.3 can nondeterministically omit complete meshes when importing
`BistroExterior.fbx` in fresh processes. The same behavior occurs with both the legacy and current
Raym0nade importer flag sets, so it is not treated as a refactor regression. Until the importer is
upgraded, patched, or replaced, verify the expected complete topology of 8,496,360 vertices and
2,832,120 faces before using a Bistro render as a visual or performance baseline. See
`docs/development-log.md` for the repeated-process evidence and current disposition.

## Dependency cleanup

Remove all project-specific native and Python packages together:

```sh
conda env remove --name raym0nade
```

This keeps Conda's shared package cache. Use `conda clean --all` separately only when no other
environment needs those cached downloads. Remove `build/` independently to discard local build
artifacts. Miniforge and the platform compiler remain separately managed system tools.

## AMD Vulkan renderer

The CPU renderer continues to use `std::thread` and remains the correctness reference. The optional
GPU backend is a headless Vulkan compute implementation selected with
`RAYM0NADE_ENABLE_VULKAN=ON` by the `gpu-debug` and `gpu-release` presets. It requires an AMD device
with Vulkan 1.2, a compute queue, buffer device address, acceleration structures, and
`VK_KHR_ray_query`; ordinary CPU builds and runs do not require Vulkan.

`VulkanPathRenderer` consumes the backend-neutral packed scene and performs tiled, deterministic
SPP accumulation with the Philox counter RNG. The shader implements the current 16-bounce path
estimator, including encoded diffuse, specular, emissive, and normal textures; alpha-tested primary,
continuation, and shadow rays; HDR environments and emissive area lights; reflection,
transmission, absorption, and nested media. It reads all G-buffer and radiance planes back into the
same `Film` representation used by the CPU renderer, so display processing and file naming share
the CPU export implementation. The smaller `VulkanPrimaryRenderer`, Ray Query intersector, probe,
and ShapeNormal benchmark remain available as focused diagnostics.

`raym0nade_gpu_render` is the general GPU rendering entry point. It accepts a checked-in console
recipe plus command-line overrides and writes the usual FXAA beauty image by default. Pass
`--all-passes` to export the full CPU-compatible Film pass set. Building the target also compiles
and validates `path_trace.comp` and copies its SPIR-V file beside the executable; normal runs do not
need a `--shader` argument.

After activating the Conda environment, configure, build, and test the optional backend with:

```sh
cmake --preset gpu-release
cmake --build --preset gpu-release
ctest --preset gpu-release
```

The configure and build commands are identical on Windows, Linux, and macOS. Run checked-in
recipes from the repository root so their relative model and output paths resolve correctly. On
Windows:

```bat
build\gpu-release\bin\raym0nade_gpu_render.exe ^
    --recipe examples\bistro_daylight_appearance_1024.txt ^
    --batch-spp 64 ^
    --output-prefix output\BistroGpuPathQuality1024
```

On Linux or macOS:

```sh
./build/gpu-release/bin/raym0nade_gpu_render \
    --recipe examples/bistro_daylight_appearance_1024.txt \
    --batch-spp 64 \
    --output-prefix output/BistroGpuPathQuality1024
```

Building is portable, but rendering still requires a device that passes the AMD Vulkan feature
gate above. See `docs/render-examples.md` for recipe overrides, output selection, and validation
examples.

Run `build/gpu-release/bin/raym0nade_gpu_probe` (`.exe` on Windows) to inspect the available device.
Add `--self-test` to build a one-triangle BLAS/TLAS and execute deterministic hardware Ray Queries;
add `--validation` during development to request the Khronos validation layer. The Conda environment
contains the Vulkan headers, loader, shader compiler, tools, and validation layers, so the
development dependencies are removed together with the environment. GPU architecture, feature
gates, packed-scene design, migration order, and performance methodology are documented in
`docs/gpu-backend.md`.

The Windows Conda layout requires the layer manifest directory to be explicit when validation is
requested. PowerShell development sessions can enable the layer and synchronization checks with:

```powershell
$env:VK_LAYER_PATH = "$env:CONDA_PREFIX\Library\bin"
$env:VK_LAYER_VALIDATE_SYNC = "true"
build\gpu-debug\bin\raym0nade_gpu_probe.exe --self-test --validation
```
