# Raym0nade

Raym0nade is an offline CPU path tracer written in C++17. It imports scenes with Assimp,
accelerates ray intersections with a BVH, and provides direct and indirect light sampling,
textured materials, HDR environment lighting, denoising, bloom, depth of field, and FXAA.

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
- `raym0nade_render_tests`: the tiny lit imported-scene render regression executable when testing
  is enabled
- `raym0nade_primary_render_tests`: the deterministic, no-file CPU primary-AOV regression
  executable when testing is enabled
- `raym0nade_console_tests`: the command-loop regression executable when testing is enabled
- `raym0nade_vulkan`: the optional experimental Vulkan backend library
- `raym0nade_gpu_probe`: the optional AMD Ray Query capability and hardware self-test executable
- `raym0nade_gpu_scene_tests`: the optional packed-scene CPU/GPU intersection comparison
- `raym0nade_gpu_primary_render_tests`: the optional CPU/GPU primary-AOV comparison
- `raym0nade_gpu_primary_benchmark`: the optional Bistro ShapeNormal comparison and GPU-only
  diagnostic utility

See `docs/architecture.md` for module boundaries, ownership rules, render data flow,
determinism guarantees, current limitations, and the staged GPU roadmap. The experimental AMD GPU
work is specified in `docs/gpu-backend.md`.

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

## GPU roadmap

The renderer currently uses `std::thread` for CPU parallelism. Experimental AMD GPU work lives on
the `dev-gpu` branch and targets a headless Vulkan compute backend with `VK_KHR_ray_query`. The
backend is optional, so ordinary CPU builds do not require a Vulkan-capable device.

The G3a/G3b vertical slices provide a shared backend-neutral render contract (`ImageExtent`,
`PinholeCamera`, `DirectionalLight`, `PrimaryRenderRequest`, `PrimaryAov`, and `LinearImage`). The
CPU implementation is a deterministic, no-file reference for `BaseColor`, `ShapeNormal`, and the
G3b `DirectDiffuse` diagnostic. The experimental `VulkanPrimaryRenderer` generates every primary
ray on the device in one two-dimensional compute dispatch using 8 x 8 workgroups and returns a
row-major linear image after one readback.

G3b adds a request-local directional light and hard opaque shadows. A miss is black; a visible hit
uses its camera-facing shape normal and evaluates
`max(baseColor, 0) * incidentRadiance * max(dot(N, L), 0) / pi`, with the light direction normalized
on the host. The same GPU invocation performs the primary query and, when needed, a terminate-on-
first-hit shadow query before writing the pixel. The render still uses one dispatch, one submit,
and one readback. `BaseColor` and `DirectDiffuse` accept only referenced opaque materials without
diffuse textures; `ShapeNormal` supports every scene accepted by the current Vulkan geometry
boundary, while alpha-cutout scenes are rejected.

G3b is a deterministic lighting diagnostic, not a complete GPU path tracer. It has no environment
or area lighting, emission, specular or metallic/roughness response, smooth normals, normal maps,
distance attenuation, random sampling, path continuation, texture sampling, accumulation, or
post-processing. The validated 4 x 4 and 13 x 9 timings are bring-up diagnostics and do not
establish a speedup; the full G3 milestone remains open.

`raym0nade_gpu_primary_benchmark` provides an explicitly limited Bistro geometry-throughput check.
It reports single-thread CPU, parallel CPU, complete GPU-call, and dispatch/readback wall-clock
timings plus a separate GPU device-timestamp duration. CPU/GPU ShapeNormal images are written
outside the timed region. Because GPU texture and candidate alpha testing are not implemented, its
benchmark-local packed-scene copy treats cutout triangles as opaque and its report records that
semantic difference. It is not a textured beauty benchmark or a substitute for the complete-
topology performance gate.

Pass `--gpu-only` for a manual GPU diagnostic that performs scene import, packing, Vulkan setup,
GPU warm-up, and measured GPU renders without executing any CPU render. This mode writes only the
GPU ShapeNormal PNG, GPU timing CSV, and an explicitly GPU-only summary. Use a dedicated output
directory so files left by an earlier CPU/GPU comparison cannot be mistaken for current output:

```sh
./build/gpu-release/bin/raym0nade_gpu_primary_benchmark --gpu-only \
    --width 3840 --height 2160 \
    --output-dir output/benchmarks/gpu-only-3840x2160
```

Windows users should invoke the corresponding `.exe`. This remains a geometry diagnostic with
benchmark-local opaque cutout fallback, not a textured or path-traced Bistro beauty render.

After activating the Conda environment, configure, build, and test the optional backend with:

```sh
cmake --preset gpu-release
cmake --build --preset gpu-release
ctest --preset gpu-release
```

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
$env:VK_LAYER_VALIDATE_SYNC = "1"
build\gpu-debug\bin\raym0nade_gpu_probe.exe --self-test --validation
```
