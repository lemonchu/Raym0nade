# Development Log

This document is the durable handoff record for the Raym0nade modernization work. Keep it
updated whenever a decision, implementation phase, validation result, or known limitation
changes.

## Project rules

- Repository content added or rewritten during the modernization must be in English.
- The CPU modernization was merged to `main`. Experimental GPU work happens on `dev-gpu` until it
  reaches its documented correctness and performance gates.
- Preserve a working CPU renderer while refactoring. GPU work starts only after the CPU
  reference path is deterministic, tested, and benchmarked.
- Project package downloads should use a direct connection rather than the local
  `127.0.0.1:7897` proxy.
- Third-party native dependencies belong to the Conda environment when practical so they can
  be removed together.
- Keep `scripts/` limited to runtime decoder modules required by the renderer. Use Conda and CMake
  commands directly instead of adding workstation-specific setup, build, or launcher wrappers.

## Current environment

- Host verified on Windows 11 x64.
- Miniforge command-line distribution, installed outside the repository
- Conda environment: `raym0nade`
- Toolchain: Visual Studio 2022 Build Tools / MSVC, CMake, and Ninja
- Main dependencies: Assimp, libpng, Python 3.12, NumPy, ImageIO, Pillow, FreeImage, and the
  optional Vulkan development stack
- Environment manifest: `environment.yml`

The environment can be removed with:

```text
conda env remove --name raym0nade
```

### AMD GPU research and Vulkan toolchain bring-up

On 2026-09-04, GPU development started from the merged CPU baseline on a new `dev-gpu` branch. A
read-only architecture audit compared HIP/HIPRT with Vulkan KHR ray tracing. The first backend is a
headless Vulkan compute implementation using `VK_KHR_ray_query`; HIPRT remains a possible A/B
candidate for a future officially supported AMD discrete GPU. This choice avoids depending on the
HIP SDK for the current unsupported `gfx1036` development APU while retaining AMD hardware ray
traversal through a ratified cross-vendor API.

The Vulkan development packages were installed directly from conda-forge without the local proxy:

- Vulkan headers and loader 1.4.357;
- shaderc 2026.3, glslang 16.5, and SPIR-V Tools 2026.3;
- Vulkan tools and Khronos validation layers 1.4.357.

The packages are declared in `environment.yml` and remain removable with the complete Conda
environment. CMake keeps the backend disabled by default and provides separate `gpu-debug` and
`gpu-release` presets. The first project-owned C++ capability probe configured and built
successfully in `gpu-debug`, then reported:

- device: integrated `AMD Radeon(TM) Graphics`, vendor/device `0x1002/0x13c0`;
- device API 1.4.315 and AMD proprietary driver dated 2026-04-17;
- compute queue family 0 and subgroup size 64;
- buffer device address, acceleration structures, Ray Query, and Ray Tracing Pipeline all enabled;
- acceleration-structure primitive limit 536,870,912; and
- successful Raym0nade AMD Ray Query capability status.

The host is a Ryzen 9 9950X with 16 CPU cores and only two integrated GPU graphics cores. It is a
valid functionality node but not a credible target for promising a large speedup over the current
CPU renderer. Performance acceptance therefore separates this integrated-GPU compatibility gate
from future RX 6000/7000/9000 discrete-GPU targets. A limited primary-AOV renderer now exists, but
no complete GPU path renderer or speedup is claimed.

The G1 hardware-intersection self-test compiles a checked-in compute shader to generated SPIR-V,
builds a one-triangle BLAS and identity TLAS, and executes deterministic hit and miss Ray Queries.
The hit returns primitive 0, distance 1, and barycentrics `(0.25, 0.5)`; the miss returns the
invalid primitive sentinel. G2 now also packs an imported scene and compares repeated Vulkan
primary-hit batches with the CPU BVH. G3a adds deterministic CPU/GPU `BaseColor` and `ShapeNormal`
images with device-side primary-ray generation. The clean validation records below close the
current Windows G1, G2, and G3a functionality slices. G3, a minimal lit GPU render, remains the next
implementation gate. See `gpu-backend.md` for the durable design and measurement plan.

## Validation history

The Bistro archive was safely extracted to `model/Bistro_v5_2/`. Model and generated output
directories remain ignored by Git.

### Legacy execution baseline

Before the modernization, the exterior daylight scene completed a `512 x 288`, 16-sample,
16-thread render. The renderer reported approximately 5.0 seconds for core rendering and 6.25
seconds total; end-to-end wall time including scene loading was approximately 18.5 seconds. The 20
outputs used the `BistroPreview512` prefix.

This run proved that the legacy asset and environment pipeline worked, but it did not capture
enough controlled camera, build, and timing metadata to serve as a pixel-exact golden image or a
rigorous before/after performance benchmark.

### Current automated validation

On 2026-09-04, both the Debug and Release build trees passed all three CTest entries on Windows 11
x64:

- `core`: numerical, geometry, distribution, texture/material, settings, selected BSDF, and BVH
  regressions;
- `render-smoke`: imports a checked-in two-triangle OBJ/MTL scene with constant diffuse and
  constant-emissive materials, renders at `8 x 8` and one sample per pixel, requires nonzero and
  equal direct-light sample counts plus visible direct diffuse lighting, and compares all 20 PNG
  outputs byte-for-byte between one and four requested worker threads; and
- `console`: verifies stream-driven command-loop behavior, including whitespace around `exit`.

The render smoke test also imports a fixture containing one valid face and one coordinate that
overflows to a non-finite value. It verifies that the damaged face is skipped and that the valid
face remains intersectable through the BVH.

The render fixture now verifies nonzero stochastic area-light transport as well as import,
G-buffer/material output, post-processing, file generation, and thread-count scheduling. It still
does not cover sky lighting, textures, transmission, cutouts, or indirect-path determinism.

The CI workflow runs the Release build and all tests on Windows, Ubuntu, and macOS, and also runs a
Windows Debug build/test job to exercise the Conda-compatible MSVC ABI configuration.

### Clean CPU and Vulkan G2 validation

On 2026-09-04, the previous generated `build/debug`, `build/gpu-debug`, and
`build/dependency-check` trees were removed before validating from clean build directories on
Windows 11 x64. No cross-platform GPU validation is implied by these results.

- CPU Debug completed all 23 Ninja build steps and all three CTest entries.
- CPU Release completed all 23 Ninja build steps and all three CTest entries.
- GPU Debug completed all 33 Ninja build steps and all five CTest entries.
- GPU Release completed all 33 Ninja build steps and all five CTest entries.
- Both Ray Query compute shaders were compiled by `glslc` for Vulkan 1.2 and accepted by
  `spirv-val`.

The clean Debug dependency database reported 170, 137, 140, and 143 dependencies for
`model.cpp`, `bvh.cpp`, `scene_data.cpp`, and `render.cpp`, respectively. This confirms that
the localized MSVC `/showIncludes` repair eliminated the earlier zero-dependency Ninja records
rather than merely reusing an old build tree.

Manual G1 runs in both Debug and Release enabled the Khronos validation layer and synchronization
validation. Each returned the expected hit and miss records with zero validation errors and zero
warnings. The latest complete manual G2 rerun used an imported packed fixture containing three
materials, six vertices, and two faces. Both repeated batches matched the CPU BVH:

| Build | Pack | Host setup | Upload | BLAS/TLAS build | Batch 1 | Batch 2 | Validation | Exit |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| Debug | 0.019 ms | 81.312 ms | 2.143 ms | 0.882 ms | 0.201 ms | 0.104 ms | layer + synchronization enabled, 0 errors / 0 warnings | 0 |
| Release | 0.003 ms | 78.993 ms | 2.285 ms | 1.040 ms | 0.197 ms | 0.135 ms | layer + synchronization enabled, 0 errors / 0 warnings | 0 |

After the timing-only test change, both `gpu-debug` and `gpu-release` were rebuilt and all five
CTest entries passed again.

These are host wall-clock bring-up diagnostics for a five-ray batch, not GPU timestamp
measurements or evidence of renderer speedup. The current device is the two-compute-unit integrated
GPU in the Ryzen 9 9950X; meaningful performance claims remain deferred to a supported AMD
discrete GPU and a complete resident rendering workload. The earlier Debug
`0x80000003` process failure did not recur in either clean configuration.

The current open G2 risks are:

- the imported comparison fixture has only two faces, so a greater-than-ten-face fixture is still
  needed to exercise CPU BVH face reordering;
- proactive Vulkan memory-budget accounting is not implemented;
- exact `tMin` and `tMax` ray-interval endpoint behavior is not covered;
- validation messages emitted during teardown are not included in the pre-destruction report;
- a fence timeout may still wait indefinitely while quiescing the queue for safe destruction; and
- GPU texture sampling and complete material shading are not implemented. G3a supports only the
  restricted constant `BaseColor` diagnostic described below, and alpha cutouts referenced by the
  packed scene are rejected rather than silently treated as opaque.

### Review hardening pass

On 2026-09-04, a clean direct CMake configuration was built and tested in both Debug and Release
on Windows 11 x64. All three registered CTest entries passed in each configuration. This pass also
confirmed that configuring the same build tree again preserves both the timestamp and SHA-256
content hash of the generated Python compatibility header; the subsequent builds performed no C++
compilation or relinking.

The corresponding fixes register the checked-in tests with CTest, reject non-finite imported face
positions before BVH construction, pass precomputed bounding-box intervals through BVH recursion,
and derive successful shutdown from the parsed console command. CI and local launcher-script entry
points were deliberately left outside this logic-focused pass.

### Script directory cleanup

On 2026-09-04, the script directory was reduced from ten source-level helpers to the two Python
modules that the current C++ runtime actually imports: `dds_to_array.py` for DDS textures and
`hdr_to_array.py` for HDR environments. Eight untracked setup, build, launch, compatibility, and
runtime-check wrappers plus generated Python bytecode were removed. Documentation, agent guidance,
and CI now invoke Conda and CMake directly, so a clean checkout has no dependency on local wrapper
scripts.

The retained decoder modules remain temporary technical debt. They must stay in place until DDS
and HDR loading is moved to native C++ or the Python bridge is packaged as a controlled runtime
resource.

The retained modules now propagate decoder exceptions to the C++ bridge, close DDS readers,
validate dimensions, channel counts, and byte-image data types, and return contiguous buffers.
Validation decoded a `2048 x 2048` four-channel Bistro DDS and the `4096 x 2048` three-channel
float32 HDR, then loaded the Bistro model successfully through the C++ bridge. Direct Debug and
Release CMake workflows both passed all three registered CTest entries without using wrapper
scripts.

### Pre-commit repository cleanup

On 2026-09-04, the generated `build/` tree and the legacy root-level `Demo.exe` and `FXAA.exe`
binaries were removed. A repository-wide scan found no batch, command, PowerShell, shell, or macOS
launcher scripts outside ignored third-party content. The `scripts/` directory contains only the
two runtime decoder modules described above, and neither module contains a workstation path.

The project was then configured from an empty build tree with the portable presets. Both
`cmake --workflow --preset debug` and `cmake --workflow --preset release` completed successfully on
Windows 11 x64, and all three CTest entries passed in each configuration. The newly generated build
tree was removed again before commit so no binary, cache, probe, or debug artifact remained in the
repository workspace.

This validation also exposed an unrelated open importer risk: two fresh Bistro processes reported
`8,441,772 / 2,813,924` and `8,496,144 / 2,832,048` vertices/faces, neither matching the previously
recorded `8,496,360 / 2,832,120` topology. The Python decoders succeeded in both runs, but the
claimed Assimp FBX topology stabilization was not holding across processes. The follow-up
investigation and decision to defer the third-party importer issue are recorded below.

### Corrected post-refactor Bistro render

The post-refactor daylight scene was rerun in two independent processes on 2026-09-04 with the
framing preserved from the `2048 x 1152` camera preset:

- Resolution: `512 x 288`
- Pixel scale: `0.0041` (four times the preset's `0.001025`, because each dimension is one quarter)
- Samples per pixel: 16
- Worker threads: 16
- Imported topology: 8,496,360 vertices and 2,832,120 faces in both processes
- Renderer-reported core render time for the recorded verification run: 4.90759 seconds
- Renderer-reported total time for the recorded verification run: 8.24438 seconds
- Direct-light samples: 1,975,269
- Outputs: 25 PNG files per run under `output/`, with prefixes `BistroCorrected512` and
  `BistroCorrected512Verify`; all 25 corresponding SHA-256 hashes match exactly

The result restores the familiar wide framing and is the current visual inspection reference, but
it is not a checked-in golden image or a controlled performance benchmark. The earlier
post-refactor `BistroRefactored512` run used the unadjusted `0.001025` scale at `512 x 288`, which
reduced both image-plane spans to one quarter of the preset and created a crop-like zoom. Its
timings remain a troubleshooting record only and must not be used as a visual baseline.

### Assimp FBX import investigation

Standalone import probes reproduced nondeterministic topology before Raym0nade loaded any material,
texture, Python decoder, project mesh, or BVH. A complete `BistroExterior.fbx` import through the
installed Assimp 5.4.3 contains 1,591 meshes, 132 materials, 8,496,360 vertices, and 2,832,120
faces. Failed runs only lost data: complete named geometry or all material-split meshes belonging to
one model disappeared together while the node count remained stable.

Each probe run used a fresh process and the same asset:

| Import configuration | Incomplete runs |
| --- | ---: |
| Current production flags and first geometry layer | 6 / 25 |
| Current flags without `PreTransformVertices` | 5 / 25 |
| No Assimp post-processing flags | 8 / 25 |
| All geometry layers and no post-processing flags | 5 / 20 |
| Legacy `32a4883` flags and default all-layer behavior | 8 / 25 |

The legacy row exactly reproduced `Triangulate | PreTransformVertices | SortByPType |
FixInfacingNormals` with `KEEP_HIERARCHY=1`. It establishes that the old Raym0nade import path also
exhibits the problem when rebuilt against the current Assimp. The historical repository did not
pin its Assimp version, so this does not prove which dependency version exhibited the behavior at
the time of the original renders.

These results rule out Raym0nade's counters, material processing, Python initialization, mesh
translation, BVH construction, `PreTransformVertices`, triangulation, and the geometry-layer option
as the direct cause. Assimp 5.4.3 also contains a known unchecked lockstep access between material
indices and face counts in `FBXConverter::ConvertMeshMultiMaterial`; upstream added a bounds check
in [commit e8651c3](https://github.com/assimp/assimp/commit/e8651c35212ab3a2b7eb97694cc7c267449bf209).
The Bistro asset was not instrumented against that source, so this is a plausible mechanism rather
than a confirmed exact cause.

The owner chose to defer the issue because the legacy import path is affected as well. The
production importer continues to select the first geometry layer, triangulate, and pre-transform
vertices, but this configuration is not considered deterministic. Bistro appearance and benchmark
runs must record and verify the complete topology before comparison. `aiProcess_SortByPType` and
`aiProcess_FixInfacingNormals` remain removed because the renderer consumes triangles after
explicit triangulation and already orients shading normals at the hit boundary.

### Higher-quality Bistro visual check

A higher-quality exterior render was completed on 2026-09-04 with these settings:

- Resolution: `1024 x 576`
- Pixel scale: `0.00205`, preserving the `2048 x 1152` preset framing
- Samples per pixel: 64
- Worker threads: 16
- Exposure: 50
- Depth-of-field post-processing: disabled
- Imported topology: 8,496,360 vertices and 2,832,120 faces
- Direct-light samples: 31,586,854
- Core render time: 58.9808 seconds
- Total render and export time: 66.627 seconds
- Outputs: 20 PNG files under `output/` with prefix `BistroQuality1024Exp50`

An initial run at the legacy exposure of 200 was visibly overexposed. Small exposure previews at
20, 50, and 100 were used to select 50 before the final run. The filtered image is substantially
cleaner than the raw estimator output, but it also shows that the current variance filter and tone
mapping remain aggressive; future image-quality work should evaluate these stages separately from
path-sampling convergence.

### Daylight appearance calibration

The apparent daylight exposure regression was traced to a preset mismatch and changed light
transport. The historical `arg_exterior` preset with exposure `200` is explicitly labelled as a
night exterior and belongs with `BistroExterior.fbx` and a `null` environment. Git contains no
complete, runnable daylight camera/exposure preset. The exterior exposure `50` record belongs to a
different camera, has no recorded environment pairing, and uses the invalid placeholder `0` SPP.
It must not be treated as a saved daylight ground truth.

Exposure application, the `0.75` luminance shoulder, hyperbolic-tangent highlight compression, and
1/2.2 display gamma are equivalent in the legacy and current implementations. For the repository's
4096 x 2048 HDR, the exact ratio between the current texel solid angle and the legacy midpoint
expression is:

```text
2 * 2048 * sin(pi / (2 * 2048)) = 3.141592346...
50 / 3.141592346...              = 15.915496...
```

The value `15.9155` corrects only that solid-angle term. It is not a complete exposure conversion:
the old sampler also conditioned samples on the shaded hemisphere without correcting its PDF, and
the current renderer adds emissive-mesh next-event estimation while an HDR is present. Both effects
vary spatially, so no global exposure produces a pixel-identical legacy image.

An intermediate 1024 x 576 render at exposure `20` completed normally with 64 SPP, 31,586,854
direct-light samples, 61.532 seconds of kernel time, and 70.0409 seconds total. The owner rejected
that result as visibly too bright. A low-resolution bracket then established exposure `14` as a
good result and requested a further exposure `12` trial. The production renderer completed that
trial with:

- Resolution: `512 x 288`
- Pixel scale: `0.0041`
- Samples per pixel: `64`
- Worker threads: `16`
- Direct-light probability: `0.7`
- Exposure: `12`
- Depth-of-field post-processing: disabled
- Output prefix: `output/BistroDaylightAppearance12Preview`
- Imported topology: 8,496,360 vertices and 2,832,120 faces
- Direct-light samples: 7,896,858
- Core render time: 14.5456 seconds
- Total render and export time: 16.409 seconds
- `Filter_FXAA` mean/median display luminance: `0.507618` / `0.490196`
- Pixels with at least one clipped display channel: `10.272%`

The runnable appearance recipes now use exposure `12` as the darker trial, with the owner-approved
exposure `14` documented as the slightly brighter alternative. `Filter_FXAA` remains the comparison
output; bloom is excluded while judging base exposure. These are scene-specific appearance
candidates, not universal legacy conversion factors.

### G3a deterministic primary-AOV vertical slice

On 2026-09-05, the first backend-neutral render boundary and device-generated primary image were
completed on `dev-gpu`. `ImageExtent`, `PinholeCamera`, `PrimaryAov`, `PrimaryRenderRequest`, and
`LinearImage` contain no Vulkan, console, file-output, sampling, or post-processing policy. The new
CPU no-file oracle renders deterministic `BaseColor` and encoded geometric `ShapeNormal` values.
It shares the exact legacy camera helper with the production CPU renderer:

```text
direction + pixelScale * ((x - width * 0.5) * right + (y - height * 0.5) * up)
```

Pixel coordinates are integers, there is no half-pixel offset, and the supplied camera basis is not
implicitly orthonormalized. Output is a row-major linear RGB image.

Packed scene format version 2 now records explicit alpha-cutout, diffuse-texture,
specular-texture, emissive-texture, and normal-texture presence bits. Texture IDs remain invalid
sentinels until device texture storage is implemented. This preserves enough source-material
capability information for a GPU backend to reject unsupported shading instead of silently using a
constant material.

The Vulkan host implementation was factored around a private `VulkanRuntime` shared by the G2
intersector and G3a renderer. It keeps device-local vertex, index, triangle-material, and material
buffers plus the BLAS/TLAS alive across operations. Acceleration-structure scratch and
instance-input buffers are freed after the build submission completes. `VulkanPrimaryRenderer`
generates the full primary-ray grid on the device in one two-dimensional compute dispatch using
8 x 8 workgroups and returns one completed `LinearImage` after one readback. It also reports an
optional GPU dispatch timestamp and handles counter wrap according to the queue family's valid
timestamp bits.

GPU `BaseColor` currently accepts only referenced opaque materials without diffuse textures.
`ShapeNormal` has the wider capability and accepts every scene that passes the current Vulkan
geometry boundary. Referenced alpha cutouts remain rejected for both AOVs.

Validation completed on Windows 11 x64:

- CPU Debug passed all 4 CTest entries.
- CPU Release passed all 4 CTest entries.
- GPU Debug passed all 7 CTest entries.
- GPU Release passed all 7 CTest entries.
- The Release primary-AOV run on `AMD Radeon(TM) Graphics` enabled the Khronos validation layer and
  synchronization validation and reported 0 errors and 0 warnings.

That Release run reported 2.861 ms for the persistent scene upload and 0.853 ms for BLAS/TLAS
construction. Its selected 4 x 4 primary diagnostics were:

| AOV invocation | Host dispatch/readback | GPU dispatch timestamp |
| --- | ---: | ---: |
| `BaseColor`, first | 0.270 ms | 0.009 ms |
| `BaseColor`, second | 0.104 ms | 0.005 ms |
| `ShapeNormal`, first | 0.093 ms | 0.006 ms |
| `ShapeNormal`, second | 0.108 ms | 0.006 ms |

This 4 x 4 test is a correctness and synchronization diagnostic, not a throughput workload or
evidence of speedup. G3a is not the complete G3 milestone: GPU lighting, random sampling, path
continuation, texture sampling, accumulation, and post-processing remain open. Performance claims
remain deferred until a representative resident workload runs on a supported AMD discrete GPU.

## Fixes completed before the modernization

- Changed the GLM submodule URL from SSH to HTTPS.
- Added the Conda environment, portable CMake presets, and an English README. Workstation-specific
  setup, build, and launch wrappers were used during bring-up but are intentionally not retained.
- Removed `Py_Finalize()` from model loading. Reinitializing NumPy-backed Python extension
  modules in the same process caused an access violation when the HDR skybox loaded after
  material textures.

## Initial audit and current disposition

The first read-only architecture, numerical-safety, performance, and build audit produced the
following findings. Their current disposition is recorded here so the original list is not
mistaken for the current implementation:

1. **Resolved:** vector finite checks now require every component to be finite.
2. **Resolved:** `ImageData` owns validated one-to-four-channel data, PNG decoding produces RGBA,
   buffer-size arithmetic is checked, and sampling cannot use an incompatible stride.
3. **Resolved:** the former image and BVH raw allocations are `std::vector` storage with RAII;
   ownership-sensitive types are non-copyable/non-movable where required.
4. **Contained, not eliminated:** `Face` still borrows material and vertex pointers. `ModelBuilder`
   computes complete capacities before filling the immutable model, and `Model`/`Bvh` cannot move,
   so the documented lifetime invariant is stable. A future packed scene should replace pointers
   with indices.
5. **Resolved:** empty and invalid random distributions return an explicit invalid sample and zero
   PDF; non-positive and non-finite weights receive no probability.
6. **Substantially resolved:** high-confidence degenerate geometry, UV/tangent, normalization,
   probability, texture-size, and medium divisions are guarded and regression-tested. Analytic
   BSDF/PDF and energy-integral testing remains open.
7. **Resolved:** model construction throws on structural/import failure, and the console publishes
   a model only after successful construction.
8. **Partially resolved:** missing normals and UVs receive defined fallbacks, non-triangle faces are
   skipped, and Assimp triangulation and vertex pre-transformation remain enabled. Bistro FBX
   topology nondeterminism is not fixed by disabling all-layer merging and remains a deferred
   third-party import risk.
9. **Resolved:** the environment is loaded before scene import, and environment and emissive-area
   lighting now coexist in the direct sampler rather than one silently replacing the other.
10. **Resolved at the API boundary:** public headers expose no Assimp, libpng, or Python types. A
    source-only `ModelBuilder` owns Assimp translation, decoder helpers are private under `src/`,
    and script/asset paths no longer depend on the process working directory. The Python decoder
    bridge now owns the GIL for each call, translates pending buffer errors, validates sizes, and
    copies potentially unaligned buffers safely. The runtime itself remains temporary technical
    debt.
11. **Resolved for the identified defect:** reflection sampling evaluates the cosine and GTR2 PDFs
    at the selected direction and includes the 50/50 proposal probability. Fresnel/refraction and
    rough-transmission terms were also repaired; broader analytic estimator validation remains.
12. **Resolved for the hot recursive path:** one row-local sample vector is reserved and reused
    across direct and recursive sampling, removing the former per-bounce vector construction and
    return-value copies. A future wavefront integrator should replace this temporary vector API
    with fixed records and queues.
13. **Resolved:** atomic row scheduling replaced the mutex queue, and stable row-derived random
    streams make the checked-in smoke fixture byte-identical across worker counts.
14. **Resolved:** the core compiles once, applications are thin targets, `Threads::Threads` is
    linked, build artifacts stay under `build/`, and shared Windows/Linux/macOS presets are present.

## Numerical and estimator status

The current primary estimator uses one Bernoulli decision per configured sample to select direct or
indirect lighting. For `0 < directLightProbability < 1`, the chosen branch is divided by its
selection probability, so low sample counts no longer lose a component through integer sample
splitting. Endpoint values remain accepted diagnostic modes: `0` disables direct lighting and `1`
disables indirect lighting, so neither endpoint is a complete beauty estimator.

The following high-confidence numerical findings are resolved:

- the RNG returns values strictly inside `(0, 1)`, and invalid distributions fail safely;
- dielectric Fresnel uses the incident angle and the actual relative IOR at entry;
- the rough BTDF includes half-vector Fresnel, Smith masking-shadowing, eta scaling, and the
  refraction Jacobian, while the macro Fresnel decision is used only as a proposal probability;
- environment next-event sampling no longer rejects valid transparent-surface directions using the
  geometric normal;
- clearcoat is white and respects material opacity;
- absorption colors are made finite and clamped per channel before exponentiation, preventing a
  color component above one from becoming path gain;
- exposure and variance finishing use double-precision scaling and finite saturation instead of
  turning overflowing highlights black; and
- failed direct or continuation proposals contribute zero instead of being retried until success,
  avoiding rejection-conditioning bias at the cost of variance.

Known deliberate biases and limitations remain:

- direct and continuation sample throughput luminance is clamped to 64 as a firefly-control
  measure;
- a factor-1 roughness floor grows monotonically along each path as regularization;
- a spatial luminance firefly clamp with a neighbor ratio threshold of 36 is applied before the
  variants named `Raw` are composed;
- paths are truncated at depth 16, and the camera uses one fixed ray per pixel rather than
  stochastic pixel-area sampling;
- the medium stack starts empty and therefore assumes the camera is outside all transmissive
  media;
- rough-dielectric lobe selection uses macro-normal Fresnel before the microfacet is sampled, so
  macro total internal reflection can remove transmission proposal support that tilted microfacets
  would still need;
- clearcoat has no dedicated GTR1 proposal and can be noisy at low roughness;
- environment importance sampling uses each selected texel's center direction, and emissive-face
  sidedness is not yet applied consistently between next-event and continuation-hit paths; and
- fixed 50/50 environment/area and cosine/GTR2 proposal choices are simple variance policies, not
  missing selection weights, but complete NEE/continuation MIS and analytic white-furnace/PDF tests
  are still absent.

## Modernization plan

### Phase 1: Safety and project boundaries - completed

- Establish a core library, thin applications, and CTest-based regression tests.
- Introduce the `raym0nade` namespace and standard public include layout.
- Replace owned raw arrays with RAII containers and make borrowed scene ownership explicit.
- Add input validation and defined failure behavior.
- Repair finite checks, RNG endpoint guarantees, distributions, image channel handling, and
  obvious divide-by-zero paths without intentionally changing the lighting model.

### Phase 2: Portable build system - completed

- Replace source globbing with explicit target source lists.
- Link imported targets for Assimp, PNG, embedded Python, GLM, and Threads.
- Keep third-party implementation details private to the core target.
- Add Windows, Linux, and macOS debug/release presets and keep artifacts inside build trees.
- Resolve Python helper scripts independently of the current working directory.

### Phase 3: Determinism and measured CPU performance - partially completed

- Seed sampling from stable render coordinates instead of worker identity.
- Replace the mutex row queue with low-overhead atomic tile or row scheduling.
- Use monotonic wall-clock timing and add repeatable smoke tests/benchmarks.
- Remove avoidable allocations from per-pixel and recursive paths where correctness can be
  preserved.

Stable row seeding, atomic row scheduling, monotonic timing, a reusable row-local sampling scratch
buffer, and a lit fixture-level cross-thread smoke test are implemented. Broader stochastic
fixtures and a controlled benchmark harness remain open.

### Later phases

- Separate scene loading, acceleration, camera, integrator, scheduler, film, post-processing,
  and image I/O into stable interfaces.
- Complete the NEE/continuation estimator design and validate BSDF/light PDFs and energy with
  analytic tests before treating appearance as a permanent contract.
- Continue the versioned packed scene and optional GPU backend beyond the completed indexed
  geometry, material-capability flags, shared primary-AOV contract, and G3a Vulkan renderer.

## Work in progress

- [x] Create and switch to the `dev` branch.
- [x] Complete the initial read-only audit.
- [x] Move application entry points into `apps/`.
- [x] Move public headers into `include/raym0nade/` and rename them to `.hpp`.
- [x] Update public APIs and namespaces after the file move.
- [x] Remove Assimp and Python implementation types from the public include tree.
- [x] Implement the modern CMake target graph and cross-platform presets.
- [x] Add numerical, geometry, distribution, texture/material, settings, BSDF, and BVH tests.
- [x] Add the tiny lit cross-thread imported-scene render smoke test.
- [x] Reproduce and scope Bistro FBX topology variation against current and legacy import flags.
- [x] Rebuild and visually inspect Bistro with resolution-correct camera scaling.
- [x] Record the available legacy and post-refactor timing observations with their limitations.
- [x] Merge the CPU refactor to `main` and create `dev-gpu` for experimental backend work.
- [x] Complete Vulkan G0/G1 and the initial imported packed-geometry G2 slice on Windows AMD.
- [x] Complete G3a deterministic CPU/GPU `BaseColor` and `ShapeNormal` primary AOVs.
- [ ] Complete G3 with GPU lighting and path integration under a CPU correctness comparison.
- [ ] Add GPU texture storage and sampling before enabling textured `BaseColor` or alpha cutouts.
- [ ] Resolve or replace nondeterministic Bistro FBX import before treating its topology as a
  deterministic regression gate.
- [ ] Add controlled before/after benchmarks, analytic estimator tests, and golden-image tolerances.

## Continuation checklist

When resuming work, run:

```text
git status --short --branch
cmake --list-presets
```

Then read this file and `README.md`. Do not assume a phase is complete until its checkbox is
updated and its validation result is recorded here.
