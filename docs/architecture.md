# Raym0nade Architecture

This document describes the architecture after the CPU renderer refactor and the first complete
Vulkan path-rendering implementation on `dev-gpu`. It is a contract for future changes, not a
description of the legacy source layout. Read it together with `development-log.md` for
chronological decisions and validation results, and with `gpu-backend.md` for the detailed AMD
GPU design.

## Architectural goals

The current design preserves a safe, portable CPU reference renderer while an optional Vulkan
Ray Query renderer is developed behind explicit capability gates. Its most important
properties are:

- one reusable C++17 core target with platform dependencies hidden behind implementation files;
- thin command-line applications instead of duplicate renderer builds;
- explicit ownership for scene, acceleration, image, and console state;
- finite-value and degenerate-input checks at numerical boundaries;
- reproducible sampling when the CPU worker count changes; and
- the same CMake presets and dependency manifest on Windows, Linux, and macOS.

The architecture remains transitional at its application and asset boundaries, but the first GPU
beauty-rendering path is complete for the current CPU estimator. `VulkanPathRenderer` consumes
immutable indexed packed data and runs the 16-bounce estimator as a per-pixel Vulkan compute
megakernel. Geometry, textures, area lights, the HDR environment, acceleration structures, and
per-tile accumulators remain device-resident across sample batches. Completed Film tiles are read
back and passed through the same host finalization and post-processing implementation as CPU
results. The older primary-AOV and arbitrary-ray paths remain as focused diagnostic interfaces.

## System overview

```text
renderer_cli
    |
    v
ConsoleApplication ---- creates/owns ----> RenderSettings
    |                                      Model
    |                                        |
    |                                        +-- Material and texture mipmaps
    |                                        +-- VertexData and Face arrays
    |                                        +-- emissive LightObject data
    |                                        +-- SkyBox importance distribution
    |                                        `-- Bvh over the Face array
    |
    `-- renderToFilm(Model, RenderSettings) -> FilmRenderResult
             |                                  |
             +-- atomic row scheduler           `-- Film buffers + RenderStats
             +-- per-row RenderContext                    |
             |   (RNG and medium stack)                    v
             `-- intersections -> BSDF/light     exportFilmToFiles
                                                   |
                                                   +-- spatial clamp and filtering
                                                   +-- shading, bloom, DOF, and FXAA
                                                   `-- PNG outputs

fxaa_cli ---- loads PNG directly into Film ---- FXAA ---- saves PNG

Backend-neutral CPU/GPU slices:

Model -----------------> CPU primary-AOV oracle(PrimaryRenderRequest) -> LinearImage
  |
  `-- packScene() -> immutable PackedSceneData -> private VulkanRuntime
                                                   |
                                                   +-- persistent geometry/light buffers
                                                   +-- BDA-paged texture texels
                                                   +-- persistent BLAS/TLAS
                                                   +-- VulkanRayQueryIntersector
                                                   |      `-- batched hit records
                                                   +-- VulkanPrimaryRenderer(PrimaryRenderRequest)
                                                   |      `-- one dispatch/readback -> LinearImage
                                                   `-- VulkanPathRenderer(RenderSettings)
                                                          +-- tiled 8 x 8 compute workgroups
                                                          +-- bounded SPP batches
                                                          +-- iterative 16-bounce transport
                                                          +-- four radiance planes + G-buffer
                                                          `-- tile readback -> Film
                                                                      |
                                                                      v
                                                            exportFilmToFiles
```

Scene construction is single-threaded. Once its constructor finishes, a `Model` is treated as an
immutable object and may be read concurrently by all render workers. Each worker writes only to
pixels in the row it has claimed.

## Build and dependency boundaries

The top-level CMake project exposes these targets:

| Target | Responsibility |
| --- | --- |
| `Raym0nade::core` | Static library containing geometry, scene loading, BVH traversal, materials, sampling, integration, film processing, and image I/O. |
| `raym0nade_renderer` | Interactive console executable, emitted as `raym0nade`; contains `renderer_cli.cpp` and the `ConsoleApplication` implementation. |
| `raym0nade_fxaa` | Standalone PNG-to-PNG FXAA utility. |
| `raym0nade_tests` | Numerical, geometry, texture, material, settings, BSDF, and BVH regression executable registered with CTest. |
| `raym0nade_render_tests` | Tiny checked-in lit-scene import/render test that compares all output passes across worker counts. |
| `raym0nade_primary_render_tests` | Deterministic no-file CPU `BaseColor`, `ShapeNormal`, and `DirectDiffuse` primary-AOV regression executable. |
| `raym0nade_console_tests` | Stream-driven command-loop regression executable, built with the console implementation and linked to the core. |
| `Raym0nade::vulkan` | Optional static library containing AMD capability discovery, a shared private Vulkan scene runtime, the Ray Query intersector, the primary-AOV renderer, and the Film-producing path renderer. |
| `raym0nade_gpu_shaders` | Optional build target that compiles the Ray Query, primary-AOV, Philox contract-test, and path-tracing compute shaders and validates their generated SPIR-V. |
| `raym0nade_gpu_probe` | Optional capability and deterministic one-triangle hardware-intersection executable. |
| `raym0nade_gpu_scene_tests` | Optional imported packed-scene CPU/GPU primary-hit comparison registered with CTest. |
| `raym0nade_gpu_primary_render_tests` | Optional deterministic CPU/GPU primary-AOV and directional-light comparison registered with CTest. |
| `raym0nade_gpu_counter_rng_tests` | Optional device-dispatch known-answer test that compares Philox integer and open-(0,1) float bits with the CPU contract. |
| `raym0nade_gpu_path_render_tests` | Optional Vulkan path-renderer regression executable covering deterministic partitioning, Film assembly, lighting, textures, cutouts, transmission, and a statistical CPU comparison. |
| `raym0nade_gpu_primary_benchmark` | Optional Bistro ShapeNormal CPU/GPU wall-clock and GPU-only diagnostic utility; generated images and reports remain under ignored `output/`. |
| `raym0nade_gpu_render` | Optional non-interactive Vulkan beauty-rendering CLI that reads one checked recipe, supports render overrides, and exports either the beauty image or all CPU-compatible Film passes. |

`Raym0nade::core` exposes the C++17 headers under `include/raym0nade/` and GLM vector types through
its build-tree interface; install/export packaging is not implemented yet.
Assimp, libpng, embedded Python, and the platform thread library are resolved through imported
CMake targets and are implementation dependencies of the core target. Assimp, libpng, and Python
types and headers do not appear in the public include tree. Assimp-specific material and mesh
translation is owned by the source-only `ModelBuilder` in `model.cpp`; decoder helper declarations
remain private under `src/`.

The console implementation is deliberately outside the core library even though its interface is
under the public include tree. This prevents an application-specific command loop from becoming a
core renderer dependency. Non-interactive frontends, including `raym0nade_gpu_render`, call scene
and rendering APIs directly and do not depend on `ConsoleApplication`.

`CMakePresets.json` uses Ninja and creates isolated `build/debug` and `build/release` trees. The
same CPU presets apply on Windows, Linux, and macOS. Separate `gpu-debug` and `gpu-release`
presets enable Vulkan when its SDK components are present; they have been validated on Windows,
not yet across all supported hosts. Conda supplies native and Python dependencies; GLM is supplied
by the pinned repository submodule. CMake does not download third-party packages.

## Module responsibilities

### Geometry (`geometry.hpp`, `geometry.cpp`)

This is the lowest renderer-specific layer. It defines vector aliases, `Ray`, `RayDifferential`,
`Box`, ray-box and ray-triangle intersection, barycentric coordinates, and tangent-space helpers.
It also owns shared numerical policy such as `kRayEpsilon`, finite checks, safe normalization,
clamped square roots, and positive-base powers.

Higher layers should use these helpers at untrusted numerical boundaries instead of duplicating
slightly different epsilon or fallback behavior.

### Components and probability (`component.hpp`, `component.cpp`)

This module defines the small data types shared by scene loading, traversal, and sampling:

- `VertexData`, `Face`, and `HitRecord` describe triangle inputs and intersection results.
- `Generator` wraps `std::mt19937` and returns floats strictly inside `(0, 1)`. Code may therefore
  use a random value as a denominator input without receiving either endpoint.
- `RandomDistribution` stores a cumulative distribution. Non-positive and non-finite weights have
  zero probability, while an unusable distribution returns `-1` and a zero PDF.
- `LightObject` groups copied emissive triangles and their face distribution.
- `SkyBox` owns an HDR environment and a luminance-times-solid-angle importance distribution.

`SkyBox::load` constructs temporary state and commits it only after validation and distribution
creation succeed. A failed environment load therefore does not leave a partially initialized
skybox.

### Acceleration (`bvh.hpp`, `bvh.cpp`)

`Bvh` builds a binary hierarchy over a mutable `std::vector<Face>`. It recursively splits the
largest centroid axis at the median with `std::nth_element`; leaves contain at most ten faces.
Interior nodes store the index of two adjacent child nodes, while leaf nodes store a range in the
reordered face array. Traversal visits the nearer child first and tightens the caller's closest-hit
distance. Each child bounding box is tested once; the resulting ray interval is passed into the
recursive visit so traversal does not repeat slab work at the child entry.

The BVH owns its node vector but does not own faces. This lifetime rule is fundamental: the face
array passed to `build` must not move, reallocate, or be destroyed while the BVH is used.

### Materials and texture data (`material.hpp`, `material.cpp`)

`ImageData` owns byte pixels for a base image and up to 32 total mip levels. Dimensions and the
channel count are private invariants exposed through read-only accessors. `setBaseLevel` accepts one
to four channels and validates dimensions, multiplication overflow, and the exact buffer size.
`generateMipmaps` uses ceil-halving and repeat-wrapped source taps for each axis, so odd and
non-square images continue to `1 x 1`. The explicit `sampleRgb` and `sampleRgba` APIs perform
wrapped bilinear/trilinear sampling; finite coordinates are reduced to the unit interval before
integer conversion, invalid coordinates return black, and invalid depth selects the base level.
One-channel data expands to opaque gray, two-channel data to gray plus alpha, three-channel data to
opaque RGB, and four-channel data remains RGBA. PNG files are decoded natively through libpng into
RGBA data with exception-safe cleanup. DDS files currently go through the private Python bridge.

`Material` owns four optional texture slots: diffuse, specular, emissive, and normal. It converts
color textures to linear space during sampling. Its public constants default to opaque, IOR 1,
roughness 0.8, metallic 0, white diffuse and transmission, and black emission; textureless
materials retain those values. Diffuse and emissive textures are multiplied by their imported
constant factors; an emissive texture without a usable imported emission receives a unit factor.
Both constant and textured emission can create emissive light objects. The current specular-map
convention takes roughness from green and metallic from blue. Assimp property interpretation and
the remaining legacy name-based transmission presets live in `ModelBuilder`, outside the public
material API.

### Scene (`model.hpp`, `model.cpp`)

`Model` is the aggregate scene owner and the boundary between asset import and rendering. A
source-only friend, `ModelBuilder`, contains all Assimp-specific work. Construction performs the
following sequence:

1. resolve and load the optional environment;
2. import and pre-transform/triangulate the model with Assimp;
3. import material constants and supported textures, resolving relative asset paths from the
   resolved model file's parent directory;
4. reserve complete vertex and face capacity;
5. create triangle, vertex-attribute, and emissive-light data; and
6. reorder the final face array and build the BVH.

Geometry or scene import failures throw and prevent the model from being exposed by the console.
Individual unsupported or failed material textures are reported and skipped so usable geometry can
still load. Imported triangles with any non-finite position are reported and skipped before BVH
construction, preventing one damaged face from invalidating otherwise usable geometry.

The FBX importer explicitly disables `AI_CONFIG_IMPORT_FBX_READ_ALL_GEOMETRY_LAYERS`. Assimp 5.4.3
was observed to omit complete meshes nondeterministically from the Bistro exterior asset, but
repeated-process probes showed that selecting the first layer does not eliminate the problem.
Removing all post-processing also reproduces it, as does the legacy importer flag set from commit
`32a4883`; the defect therefore precedes Raym0nade's material, Python, mesh-conversion, and BVH
work. The first-layer setting remains because the renderer does not need duplicate geometry layers,
but it is not considered a stabilization fix. Resolution is deferred; representative
repeated-process checks and the expected topology manifest are required when changing Assimp
versions or importer flags. Vertex pre-transformation remains required because the current packed
CPU faces do not retain or evaluate the source node hierarchy.

`Model::intersect` and `Model::occluded` apply alpha-cutout traversal on top of the BVH, with a
maximum of 32 transparent layers. Hit material evaluation interpolates UVs and normals, estimates a
texture footprint from ray differentials, and applies a guarded tangent-space normal map.

### Backend-neutral primary rendering (`render_contract.hpp`, `render.hpp`)

`ImageExtent`, `PinholeCamera`, `DirectionalLight`, `PrimaryAov`, `PrimaryRenderRequest`, and
`LinearImage` form the first renderer-facing contract shared by the CPU and Vulkan implementations.
The request is independent of file naming, post-processing, threading, random sampling, and
backend SDK types. `DirectionalLight` stores a direction from the shaded point toward an infinitely
distant light plus its linear incident radiance; validation requires a finite nonzero direction and
finite nonnegative radiance. `LinearImage` owns row-major linear RGB values indexed as
`y * width + x`.

The camera helper deliberately preserves the existing renderer's integer-pixel convention. For
integer coordinates `(x, y)`, the normalized primary direction is derived from:

```text
direction + pixelScale * ((x - width * 0.5) * right + (y - height * 0.5) * up)
```

There is no half-pixel offset, and the supplied `right` and `up` vectors are not silently
orthonormalized. The legacy CPU renderer and CPU primary-AOV oracle call this helper; the Vulkan
shader implements the same contract and is checked against the oracle, so camera changes cannot
drift between those paths unnoticed.

`renderPrimaryAovCpu` is the deterministic, no-file correctness oracle. It traces exactly one
primary ray per pixel and returns evaluated `BaseColor`, encoded geometric `ShapeNormal` value
`(normal + 1) * 0.5`, or G3b `DirectDiffuse`; misses are black. `DirectDiffuse` orients the shape
normal toward the camera, normalizes the request direction, traces a conditional hard-shadow ray,
and evaluates:

```text
max(baseColor, 0) * incidentRadiance * max(dot(N, L), 0) / pi
```

It performs no random sampling, post-processing, or image export, and repeated calls with the same
model and request are pixel-identical. The two-argument overload remains a single-thread oracle.
`CpuPrimaryRenderOptions` is a CPU-only scheduling argument for the three-argument overload: a
positive count selects that many workers, zero selects hardware concurrency, and the resolved count
is clamped to the image height. Every worker claims complete rows from a relaxed atomic counter, so
one-worker and multi-worker images are component-exact for the current primary AOVs.

### Packed scene and Vulkan execution (`scene_data.hpp`, `gpu/`)

`Model::packScene` derives a backend-neutral value representation without changing the immutable
CPU model. `PackedSceneData` owns 32-byte aligned vertex records, tightly packed 32-bit triangle
indices, one material ID per triangle, fixed-size material records, deduplicated texture
descriptors, complete mip descriptors, encoded RGBA8 texel words, indexed area-light records, and
HDR environment records. Compile-time ABI checks and runtime finite-value, count, flag,
reserved-field, range, and index validation protect the host/shader boundary.

Packed format version 4 gives every present diffuse, specular, emissive, or normal texture a
validated ID. Normalized source paths are compared without case folding. Each area light stores its
center and power plus one contiguous range of triangle records; each triangle stores three packed
vertex IDs, its emissive material ID, geometric area, normalized probability, and CDF. The
environment stores linear RGB radiance, exact per-row texel solid angles, a row distribution, and a
conditional per-row texel CDF. Packing and validation reproduce the CPU importance weights from
luminance multiplied by exact solid angle, verify emissive geometry against packed positions, and
reject malformed ABI metadata or distributions.

The private `VulkanRuntime` is shared implementation infrastructure for the arbitrary-ray,
primary-AOV, and path-rendering APIs. Each owning backend object selects a compatible AMD compute
device, uploads immutable geometry, material, texture metadata, area-light, and environment arrays
to persistent device-local buffers, and builds a persistent BLAS plus identity TLAS. A reusable
staging allocation transfers large arrays in chunks no larger than 16 MiB. Acceleration-structure
scratch and instance-input buffers are released once the build submission completes. Vulkan
handles and SDK types remain below the public include tree.

Encoded texture texels are not constrained to one storage-buffer descriptor. The runtime divides
the logical `uint32_t` texel array into power-of-two pages, normally no larger than 256 MiB and
always bounded by `maxStorageBufferRange`. A compact descriptor-visible table stores the page
size, shift, count, total texel count, and each page's 64-bit Vulkan buffer device address. Shader
sampling resolves a logical texel ID to a page and offset, then reads that page through
`GL_EXT_buffer_reference`. The runtime validates address alignment, page ranges, the table size,
the `uint32_t` scene limit, and `maxMemoryAllocationCount` before publication.

`VulkanRayQueryIntersector` reuses that runtime for complete ray batches and returns hit flag,
primitive ID, parametric distance, and barycentrics through Vulkan-free public records.
Its compact ray-batch API deliberately rejects cutout scenes because it has no texture-footprint or
material-sampling contract.
`VulkanPrimaryRenderer` instead generates every camera ray on the device and executes one
two-dimensional compute dispatch using 8 x 8 workgroups. For `DirectDiffuse`, the same shader
invocation reuses one `rayQueryEXT` variable for the primary query and then, only for a front-lit
hit, for a shadow query. The host normalizes the light direction and
supplies it with the incident radiance through a 112-byte push-constant block. The host also stores
`std::nextafter(kRayEpsilon, +infinity)` in the otherwise padded `cameraDirection.w`; both primary
and shadow queries use that value as `tMin`. This maps Vulkan's closed lower bound to the CPU
intersector's strict `t > kRayEpsilon` rule without changing the push-constant size. G3b did not
change the packed-scene ABI or descriptor bindings. Each render still uses one dispatch, one
submission, and one completed row-major `LinearImage` readback.

Repeated renders reuse the packed scene buffers and acceleration structures. GPU timestamp queries
measure the compute dispatch when supported, including correct wrap handling for the queue family's
reported valid-bit count; host timing separately includes submission, waiting, and readback.

`shaders/include/packed_scene.glsl` is the shared GPU layout, texture, and traversal implementation.
The primary shader samples diffuse textures for `BaseColor` and `DirectDiffuse` with the CPU's V
flip and repeat wrap, bilinear/trilinear filtering over encoded mip bytes, then RGB
`pow(value, 2.2)` decoding. Diffuse footprints come from primary-ray differentials. Alpha-cutout
traversal samples base-mip alpha and rejects values below `1e-4` in the Ray Query candidate loop for
both primary and shadow rays.

Vulkan does not guarantee Ray Query candidate order. The GPU loop therefore has no traversal-order
counter: it rejects every transparent candidate and confirms every acceptable triangle, allowing
the implementation to choose the nearest confirmed hit. This preserves normal cutout semantics and
avoids false misses caused by an arbitrary 32-candidate limit. It intentionally differs only for a
pathological ray crossing more than 32 distance-ordered transparent layers: the CPU's defensive cap
returns a miss, while the GPU continues traversal.

The primary renderer remains a deterministic diagnostic. `VulkanPathRenderer` is the complete
beauty-rendering path for the current CPU estimator. It accepts the existing `RenderSettings`,
uses one invocation per pixel within each tile, and runs an iterative transport loop capped at 16
bounces. It reproduces the direct/indirect Bernoulli estimator, sixteen-replicate rule for a
transparent first hit, roughness continuation rule, throughput clamp, reflection and transmission
branches, nested-medium IOR and absorption state, emissive endpoints, and environment/area-light
next-event estimation.

The path shader evaluates diffuse, specular, emissive, and normal textures using the packed
material convention; diffuse alpha drives candidate-confirm cutout traversal. It reconstructs the
first-hit G-buffer and accumulates direct diffuse, direct specular, indirect diffuse, and indirect
specular radiance. Each radiance record carries its RGB sum and scalar second moment so host
readback can construct the same `RadianceData` representation as the CPU renderer.

SPP is divided into configurable batches of at most 64 samples. The output allocation remains
resident for every batch of one tile. Each batch is submitted and fence-waited separately, with a
shader-write to shader-read/write dependency before accumulation resumes; this bounds an individual
Windows submission rather than placing all SPP in one watchdog-sensitive command buffer. After the
last batch, a separate transfer submission copies the tile to host memory. Tile size and batch size
change scheduling only: Philox addresses remain keyed by seed, image pixel, sample, bounce,
replicate, and dimension.

The host assembles readback tiles into a `Film`, validates shader-produced IDs and Boolean fields,
and calls `finalizeRadianceData` for each of the four accumulated planes. The returned
`VulkanPathRenderResult` is therefore compatible with `exportFilmToFiles` and with the same
spatial clamp, filtering, shading, bloom, depth of field, FXAA, and PNG code used for CPU results.

### Sampling (`sampling.hpp`, `sampling.cpp`)

`Bsdf` evaluates and samples the current diffuse, microfacet reflection, clearcoat, and rough
transmission model. Dielectric Fresnel uses the incident angle and relative IOR, and the rough BTDF
includes its microfacet distribution, masking-shadowing term, eta scaling, and refraction Jacobian.
Reflection sampling uses a 50/50 cosine/GTR2 proposal and evaluates both proposal PDFs at the
selected direction. Direct-light sampling chooses between the environment and emissive triangle
categories, using a fixed 50/50 choice when both exist, then performs visibility testing through
`Model::occluded`. Those category and reflection probabilities are included in the corresponding
sampling densities; they are variance choices rather than unaccounted estimator weights. Invalid
sampling attempts are not appended, but the averaging denominator remains the requested sample
count, so they are mathematically zero-contribution trials rather than a distribution conditioned
by retrying until success.

`LightSample` keeps estimator terms explicit:

- `throughput` is the BSDF/path response;
- `radiance` is the incident radiance including the inverse sampling density used by the sampler;
- `weight` is the averaging weight assigned by the calling estimator.

Keep these meanings stable. Reusing these fields as anonymous color accumulators makes estimator
reviews and GPU translation unnecessarily difficult.

### Integration (`render.hpp`, `render.cpp`)

`RenderSettings` is the validated render request: camera frame, image dimensions, sampling mixture,
thread count, deterministic seed, and output prefix. A zero thread count means automatic hardware
concurrency; the resolved count is clamped to the image height.

The legacy `pixelScale` setting is an image-plane/angular scale **per pixel**, not a
resolution-independent field of view. Preserving the same framing while scaling both image
dimensions by a factor `s` requires dividing `pixelScale` by `s`. For example, the Bistro exterior
preset uses `0.001025` at `2048 x 1152`, so the equivalent `512 x 288` preview uses `0.0041`.
Replacing this coupling with an explicit vertical or horizontal field-of-view camera parameter is
a future public-API improvement.

`renderToFilm` is the complete no-file CPU renderer. It traces one camera ray per pixel, stores
the first visible surface in the G-buffer, and treats each configured sample as a Bernoulli mixture
of first-hit direct lighting and a recursive indirect path. For a probability strictly between zero
and one, each selected branch is divided by its selection probability; this avoids the low-sample
integer-rounding loss of either component. A value of zero or one deliberately disables one branch
and is therefore useful for component diagnostics, not a complete beauty estimate. Transparent
first hits run sixteen indirect replicates and average them within the selected indirect sample.
The integrator supports reflection, transmission, absorption, and a nested-medium stack. Recursive
paths are capped at 16 bounces and use roughness-based continuation/direct-light decisions.

The returned `FilmRenderResult` contains the untouched integration buffers and render-only
statistics. `exportFilmToFiles` consumes a film value, applies the spatial high-variance clamp,
produces the variants named `Raw`, filters the four radiance components, and exports the existing
display, bloom, depth-of-field, and FXAA variants. Consequently, `Raw` still means unfiltered after
the spatial clamp, not untouched estimator output. `renderToFiles` is the compatibility wrapper
that retains progress messages, total timing, fixed pass naming, and PNG output.

### Film and post-processing (`image.hpp`, `image.cpp`)

`Film` owns all full-frame storage through `std::vector`:

- a `HitInfo` G-buffer;
- direct/indirect diffuse and specular `RadianceData` buffers; and
- the reusable display/output pixel buffer.

It validates dimensions and buffer sizes before processing. Its operations include spatial
outlier clamping, geometry-aware radiance filtering, composition of diagnostic or beauty passes,
bloom, depth-of-field blur, gamma/tone mapping, FXAA, and native PNG load/save. Filtering the four
radiance buffers uses four independent threads; each thread mutates a different buffer and reads
the shared G-buffer.

`finalizeRadianceData` centralizes the exposure and variance finishing rule for one
`RadianceData` sample accumulator. Both the CPU integrator and Vulkan Film readback use this
public host-side helper, avoiding separate exposure, second-moment, overflow, or non-finite
behavior.

`saveLinearImagePng` is the narrow display/export path for a validated `LinearImage`. It applies an
optional nonnegative linear scale, then the same highlight shoulder and gamma transform as `Film`,
and allocates only the final three-byte-per-pixel RGB buffer. This avoids constructing the full
G-buffer and radiance-buffer set merely to inspect a CPU/GPU diagnostic image.

### Console and applications (`console.hpp`, `console.cpp`, `apps/`)

`ConsoleApplication` owns named models as `std::unique_ptr<Model>` and named settings by value. It
injects input/output streams, which keeps parsing separate from process-global standard streams and
makes command behavior testable. A newly constructed model is inserted only after its constructor
succeeds.

The interactive renderer and FXAA `main` functions contain no renderer policy. The GPU rendering
CLI is intentionally non-interactive: it reads one single-model/single-settings console recipe,
applies validated resolution, SPP, seed, exposure, probability, output, tile, and batch overrides,
constructs one packed scene and `VulkanPathRenderer`, and exports a beauty image by default or all
Film passes on request. CMake copies the generated path-tracing SPIR-V beside the executable; an
explicit shader override and a build-tree fallback remain available.

## Ownership and lifetime contracts

The scene intentionally uses some borrowed pointers for CPU traversal performance. Their validity
depends on the following invariants:

| Object | Owns | Borrows / exposes | Required invariant |
| --- | --- | --- | --- |
| `ConsoleApplication` | `Model` instances and `RenderSettings` maps | Streams supplied to its constructor | Streams outlive the application. |
| `Model` | materials, vertex data, faces, lights, skybox, and BVH object | Read-only access to lights and skybox | No scene mutation after construction. |
| `Face` | three positions | Three `VertexData` pointers and one `Material` pointer | Referenced model arrays must remain at stable addresses. |
| `Bvh` | node vector | Pointer to `Model::faces_` storage | Face storage is not resized or moved after `build`. |
| `LightObject` | copied `Face` values and its distribution | The copied faces retain vertex/material pointers | The parent `Model` remains alive. |
| `HitRecord` | distance range | Pointer to a face in the model's reordered face array | Consume only while that model remains alive. |
| `HitInfo` | A value snapshot of evaluated surface data | Nothing | Safe to store in the film independently of a face pointer. |
| `Film` | Every image/G-buffer/radiance allocation | Public references obtained by callers | Do not resize buffers while rendering or filtering. |
| `LinearImage` | Row-major linear RGB pixels and their extent | Nothing | Validate extent and pixel count before publication. |
| `PackedSceneData` | Indexed vertices, triangles, material IDs, materials, texture descriptors/mips/texels, area-light triangles, and HDR radiance/distributions | Nothing | Validate before publication; treat as immutable while any backend upload is derived from it. |
| `FilmRenderResult` | One complete no-file CPU `Film` and its render-only statistics | Nothing | Export by moving or copying the result film; export processing intentionally mutates its private value. |
| Private `VulkanRuntime` | Vulkan instance/device, queue synchronization, persistent device-local scene buffers, BLAS, and TLAS | Packed scene data only during synchronous construction | Serialize command-buffer use; submitted work must quiesce before owned resources are destroyed. |
| `VulkanRayQueryIntersector` | One `VulkanRuntime`, batched-ray buffers, pipeline, and descriptor state | Nothing after construction | Do not call one intersector concurrently. |
| `VulkanPrimaryRenderer` | One `VulkanRuntime`, reusable output/readback buffers, timestamp query pool, pipeline, and descriptor state | Nothing after construction | Do not call one renderer concurrently; the packed scene's supported-feature restrictions are checked per requested AOV. |
| `VulkanPathRenderer` | One `VulkanRuntime`, reusable tile output/readback buffers, timestamp state, path pipeline, and descriptor state | Nothing after construction | Do not call one renderer concurrently; scene buffers persist across renders, while each render serializes tile batches and returns an independent Film. |

Before mesh conversion, `Model` reserves the total Assimp-reported vertex and face counts. Material
storage is sized once. These steps keep pointers embedded in `Face` stable. Both `Model` and `Bvh`
are non-copyable and non-movable so an accidental container operation cannot silently violate the
BVH's borrowed face pointer. Any future mutable-scene feature must replace these pointer contracts
with stable indices or rebuild all dependent structures after mutation.

## Render data flow

For each pixel, the CPU path performs these stages:

1. Build and normalize a camera direction and its direction differentials.
2. Intersect the immutable model through the BVH and alpha-cutout layer handling.
3. Store the first-hit `HitInfo` in the film G-buffer, including material and normal-map results.
4. For every configured sample, select the direct or indirect estimator stochastically and apply
   inverse selection-probability compensation.
5. For direct samples, sample an environment or emissive triangle, test visibility, evaluate the
   BSDF, and accumulate diffuse/specular radiance and variance.
6. For indirect samples, sample reflection or transmission, update the medium stack, recursively
   trace until a light/environment contribution, termination decision, invalid value, or depth 16.
7. Store radiance components in the pixel's four film buffers.
8. After all rows join, derive and save the requested diagnostic and beauty outputs from the same
   G-buffer and radiance data.

The Vulkan path follows the same estimator stages within one compute invocation per tile pixel.
It replaces CPU recursion with a bounded loop and replaces the row-local generator with addressed
Philox samples. The first SPP batch writes the G-buffer and initializes four device accumulators;
later batches resume those accumulators. Only a completed tile crosses back to the host, where its
G-buffer and radiance planes are assembled into the common `Film`.

The renderer sanitizes invalid directions, PDFs, radiance, and throughput at several boundaries.
Invalid samples contribute zero rather than contaminating an entire image with NaNs. This is a
safety policy, not a substitute for tests proving that ordinary samples are valid.

## Threading and determinism

Rendering uses an atomic row counter. A worker claims the next row with a relaxed atomic increment,
renders every pixel in that row, and then claims another row. There is no mutex-protected work
queue and no shared random generator.

Each row creates one `RenderContext` containing:

- a `Generator` seeded by a stable hash of `RenderSettings::seed` and the row index;
- a path-local `MediumStack`; and
- a reserved `LightSample` scratch vector reused by direct and recursive sampling; and
- a local direct-light sample counter.

Because a row's seed and its left-to-right sequence do not depend on which worker claims it,
changing `threadCount` does not change sampling order within any pixel. Workers write disjoint film
rows, and only the aggregate sample counter is atomic. The checked-in render smoke test verifies
byte-identical output for all 20 PNG passes with one and four requested workers on its tiny scene.
The fixture contains one constant-diffuse triangle and one constant-emissive area-light triangle;
it also requires nonzero, equal direct-light sample counts and visible direct diffuse lighting.
This exercises deterministic stochastic direct-light transport in addition to import, G-buffer
output, post-processing, and scheduling. Broader output is expected by design to be bit-identical
across worker counts, although sky, textured, transmission, and indirect-path fixtures are still
needed. Exact equality is not promised across compilers, architectures, math-library versions, or
build configurations.

Changes that add stochastic work must preserve the mapping from stable render coordinates to
random streams. The Vulkan Philox contract uses
`(seed, pixel, sample, bounce, replicate, dimension)` and is independent of tile size, SPP batch
size, dispatch order, or worker scheduling. Its integer output and open-(0,1) float conversion have
a shared CPU/GLSL known-answer contract. Never seed from a worker index, clock, or atomic work
order. Floating-point path results are repeatable for a fixed device, driver, shader binary,
settings, and seed; bit identity is not promised across different GPU implementations.

The four-thread film filter is deterministic because each worker owns one complete radiance buffer.
Other post-processing stages currently run after render workers have joined and are mostly
single-threaded.

The CPU primary-AOV overload uses the same disjoint-row scheduling rule but has no random state.
Its `ShapeNormal` branch evaluates geometric and interpolated normals directly after the hit and
does not calculate ray differentials or sample unrelated material channels. Alpha-cutout handling
still occurs inside `Model::intersect` before the AOV-specific hit evaluation.

## Failure and validation boundaries

- Invalid render dimensions, sample counts, camera data, probabilities, or output paths throw
  before worker creation.
- Model import and structurally invalid meshes throw during construction.
- Missing or failed individual textures are non-fatal and leave the corresponding slot empty.
- Empty probability distributions fail with an explicit sentinel instead of division by zero.
- Intersection functions are `noexcept`; malformed or degenerate inputs produce a miss.
- PNG and Python-decoded buffers are checked against dimensions and channel counts before copying.
- Output directory creation and image writes propagate filesystem/libpng exceptions to the console.
- Vulkan constructors reject unsupported devices, incompatible scene ABI, invalid buffer/address
  ranges, and unsupported allocation counts before exposing a renderer.
- GPU Film readback rejects invalid shader-produced material IDs and Boolean fields, then applies
  the common non-finite and overflow sanitization while finalizing radiance.

Library callers should treat these exceptions as operation failures. Do not continue with a model
whose construction failed, and do not suppress render validation errors inside worker threads.

## Current limitations and debt

These are known boundaries of the current implementation, not accidental omissions from this
document:

1. The CPU is the only portable renderer backend. The experimental Vulkan module is restricted to
   compatible AMD devices, but its `VulkanPathRenderer` implements the complete current
   16-bounce CPU estimator and returns a CPU-compatible `Film`. Post-processing and PNG export
   remain on the CPU; there is no GPU post-processing kernel or common polymorphic backend
   interface.
2. The complete CPU integrator now returns a no-file `FilmRenderResult`, but its `RenderSettings`
   contract still includes an output prefix and the API does not yet expose cancellation,
   selectable AOV export policy, or a polymorphic backend interface.
3. CPU traversal still consumes pointer-rich array-of-structures data: faces duplicate positions,
   light objects copy faces, and the BVH points into host memory. `Model::packScene` derives a
   compact indexed version-4 upload representation with textures, area lights, and HDR environment
   data, and the Vulkan path consumes all of those arrays. The CPU renderer does not yet consume
   the packed representation directly, so CPU/GPU parity is maintained across two implementations
   rather than one shared shading kernel.
4. BVH construction uses a median split rather than SAH and has no instancing, refit, parallel
   build, wide-node layout, or stackless/device traversal form. Assimp pre-transforms vertices, so
   source hierarchy and instances are flattened during import. The Bistro FBX topology has also
   varied across fresh processes despite disabling all-layer merging, so deterministic import is
   not yet guaranteed for that asset.
5. DDS textures and HDR environments still depend on a persistent embedded Python interpreter and
   repository decoder scripts. The bridge validates and copies Python buffers and acquires the GIL
   for each call, but concurrent full-model loading has no explicit contract or regression test.
   Native decoders should replace this bridge before a standalone renderer SDK or GPU asset
   pipeline is considered complete.
6. Material interpretation supports only the current four texture slots, the first texture of each
   slot, binary opaque/transmissive handling for imported opacity, and several legacy name-based
   transmission presets. Surface specular remains fixed at 0.04 and clearcoat parameters are not
   exposed as material properties. Mip generation averages encoded bytes rather than linear-light
   color, and specular maps currently sample only their base level. This is not a general
   glTF/Disney material implementation.
7. The light and BSDF estimators remain a reference-in-progress. Fixed 50/50 environment/area and
   cosine/GTR2 proposal mixtures include their selection probabilities and are primarily variance
   choices, not missing estimator weights. However, rough-dielectric lobe selection uses Fresnel at
   the macro surface normal before sampling a microfacet: macro total internal reflection can assign
   zero transmission probability even when tilted microfacets have nonzero BTDF support, biasing
   those cases. Clearcoat has no dedicated GTR1 proposal, so sharp clearcoat can have high variance.
   Environment importance sampling selects a texel by luminance and solid angle but samples only its
   center direction. Area NEE accepts only front-facing emissive triangles, while continuation-hit
   emission has no equivalent facing check; the intended emission-sidedness contract still needs to
   be defined and tested. There is no complete MIS strategy joining next-event and continuation
   sampling. The current row-local vector scratch avoids per-bounce allocation, but a wavefront
   integrator still needs fixed records and queues. Analytic white-furnace, PDF-normalization, and
   estimator-integral tests are still required.
8. Preview-oriented firefly controls deliberately bias the current result: direct and continuation
   sample throughput luminance is clamped to 64, the minimum roughness grows monotonically along a
   path with a regularization factor of 1, and a spatial luminance firefly clamp with a neighbor
   ratio threshold of 36 is applied before even the outputs named `Raw` are composed. The maximum
   depth of 16 also truncates remaining
   transport. In addition, the camera traces one fixed ray through each pixel, so samples do not
   integrate pixel area and FXAA is the only current anti-aliasing stage. A selectable reference
   mode should disable the throughput clamp and path regularization and expose pre-spatial-clamp
   radiance before CPU output is treated as an unbiased oracle.
9. `directLightProbability` is correctly compensated only when it is strictly between zero and
   one. The accepted endpoint values deliberately remove the indirect or direct component and must
   be treated as diagnostic/AOV-like modes rather than complete beauty renders.
10. Camera validation requires finite, nonzero, linearly independent axes but does not require an
    orthonormal frame. The legacy per-pixel `pixelScale` also couples framing to resolution until a
    field-of-view API replaces it.
11. The medium stack always starts empty and assumes the camera begins outside transmissive media.
    A camera placed inside a dielectric receives incorrect initial IOR and absorption state.
12. A fixed world-space ray epsilon may fail on scenes with very different scales. CPU
    alpha-cutout traversal also has a 32-layer defensive cap; Vulkan candidate traversal has no
    ordering-dependent cap and therefore differs only for unusually deep transparent stacks.
13. Film memory is several full-resolution buffers, and bloom/depth-of-field/FXAA allocate
    additional full-frame temporaries. Post-processing and writing many PNG variants can dominate
    small renders.
14. Current tests cover numerical helpers, distributions, intersections, tangent bases, BVH
    behavior, texture and material invariants, selected dielectric behavior, render-setting
    validation, a tiny lit imported-scene direct-light render, and deterministic CPU/GPU primary
    `BaseColor`, `ShapeNormal`, and directional `DirectDiffuse` AOVs. Vulkan fixtures additionally
    exercise Philox known answers, path-render tile and SPP-batch invariance, Film assembly,
    environment misses, emissive hits, area-light transport, all four material texture slots,
    cross-page texture reads, alpha cutouts, two-boundary transmission, mixed direct/indirect
    transport, and a high-SPP statistical CPU comparison. The project still needs broader loader
    fixtures, analytic BSDF/PDF and energy tests, stored beauty-image goldens with documented
    tolerances, sanitizer runs, and measured cross-platform performance gates.
15. Vulkan development and runtime coverage remain limited to the current Windows AMD integrated
    GPU; macOS and Linux GPU builds and discrete AMD hardware remain unverified. The imported
    intersection fixture has two faces and does not exercise the CPU BVH's greater-than-ten-face
    reorder path. Proactive heap-budget enforcement, teardown-time validation capture, G2
    arbitrary-`tMax` endpoint coverage, and bounded fence-timeout recovery remain open.
    Texture texels now cross the single-buffer limit through buffer-device-address pages, but the
    full decoded RGBA8 texel array still resides in host packed data during setup and in device
    pages afterward. There is no native compressed image storage, streaming, eviction, or
    out-of-core rendering, and the logical texel array remains limited to `uint32_t` indexing.
    Extremely small float probabilities can collapse adjacent CDF entries into plateaus; selection
    code tolerates zero-width bins, but their precision loss remains a packed-format limitation.
    Making the complete BLAS candidate geometry whenever any cutout is present may also reduce
    traversal performance.
16. The per-pixel megakernel has no wavefront queues or compaction. Per-batch submission and fence
    waiting bounds watchdog exposure but adds synchronization overhead, and each completed tile is
    copied to the host before the next tile. The current integrated AMD functionality device has
    not demonstrated a speedup over the CPU reference; performance work requires profiling on a
    supported discrete AMD GPU rather than assuming feature parity implies acceleration.

## GPU-backend roadmap

GPU work should preserve the CPU renderer as a correctness oracle and proceed through explicit
compatibility gates.

### Stage 0: Lock the CPU reference

- Extend the existing dielectric spot checks with analytic tests for BSDF energy/PDF consistency,
  light estimators, refraction, medium nesting, texture filtering, and cutout traversal.
- Promote the checked-in tiny-scene smoke test from cross-thread byte comparison to stored golden
  passes with documented tolerances and add more loader/material fixtures.
- Add a selectable reference mode that disables throughput clamping and path roughness
  regularization, and make the diagnostic semantics of mixture-probability endpoints explicit.
- Benchmark scene load, BVH build, trace time, post-processing time, peak memory, and rays/s.
- Resolve or explicitly bound remaining estimator bias before treating CPU images as the permanent
  visual contract.

Gate: Debug and Release tests pass on Windows, Linux, and macOS; fixed-seed output matches across
CPU thread counts; benchmark inputs and measurements are reproducible.

### Stage 1: Introduce stable renderer contracts

- Split `renderToFiles` into scene loading, `RenderRequest`, renderer execution, film result, and
  export policy.
- Add a small backend interface with an explicit capability query and retain a
  `CpuRenderBackend` implementation.
- Make camera, integrator settings, output layers, cancellation, and statistics explicit request
  data rather than console behavior.
- Keep all file I/O and interactive parsing above the backend boundary.

Status: substantially implemented. The backend-neutral extent, pinhole-camera, directional-light,
primary-AOV request, and linear-image contracts support deterministic no-file CPU and Vulkan
primary rendering. Complete CPU integration returns `FilmRenderResult`; `VulkanPathRenderer`
accepts the same validated `RenderSettings` and returns a Film-bearing result; Film export is a
separate consuming operation for both. Cancellation, explicit output-layer policy, and a common
polymorphic backend interface remain open.

Gate: the existing CLI and golden CPU images use the new interface without visual changes.

### Stage 2: Build a device-ready scene representation

- Replace intra-scene pointers with 32-bit or 64-bit indices and immutable spans/views.
- Pack positions, attributes, material parameters, texture descriptors, lights, and flattened BVH
  nodes into versioned arrays with checked size/alignment rules.
- Preserve an import-side `SceneAsset` and derive separate CPU/GPU-ready `SceneView` data so Assimp
  never enters a render kernel boundary.
- Replace Python DDS/HDR decoding with native asset decoding and define a consistent linear-color,
  alpha, mipmap, and environment convention.

Status: complete for the Vulkan path boundary. Packed format version 4 supplies validated indexed
geometry, triangle material IDs, fixed material records, deduplicated texture IDs, complete mip
descriptors, encoded RGBA8 texel storage, indexed area lights, and HDR radiance with a hierarchical
importance distribution. Vulkan uploads and consumes all of these arrays. Buffer-device-address
pages allow the encoded texel array to exceed one storage-buffer range. A packed-scene CPU shading
path, native compressed-image storage, streaming, and proactive heap-budget enforcement remain
open.

Gate: packed-scene CPU traversal and shading match the Stage 0 reference before any GPU kernels are
accepted.

### Stage 3: Convert the integrator to device-ready iterative execution

- Replace recursive `std::vector<LightSample>` returns with fixed records and an iterative path
  state machine. Validate that transformation on the CPU before porting it.
- Use a counter-based random sequence keyed by pixel, sample, bounce, and dimension. This removes
  dependence on queue compaction and GPU scheduling order.
- Define fixed capacities, overflow behavior, termination flags, and per-path medium state.
- Start the GPU implementation with one invocation per pixel, small SPP dispatch batches, and an
  iterative megakernel. Introduce compacted wavefront queues only if profiling demonstrates that
  divergence or register pressure justifies their synchronization and storage cost.

Gate: the iterative GPU result agrees statistically with the recursive CPU reference, and random
streams remain deterministic when execution partitioning changes.

Status: implemented through the GPU path while the CPU reference intentionally remains recursive.
A shared Philox4x32-10 CPU/GLSL contract maps seed, pixel, sample, bounce, replicate, and dimension
to endpoint-safe values independently of scheduling. The Vulkan shader uses fixed 17-entry medium
state, an iterative loop through depth 16, addressed direct-light samples, and per-pixel Film
accumulators. Tile and SPP-batch partitioning do not change random addresses. Replacing CPU
recursion with the same state machine is no longer required for GPU execution, but remains useful
if future work needs instruction-level CPU/GPU comparison.

### Stage 4: Add the optional Vulkan Ray Query backend

- Use headless Vulkan compute with `VK_KHR_ray_query`, initially gated to AMD devices. Keep Vulkan
  types out of public scene and renderer headers. A full Ray Tracing Pipeline and HIPRT remain
  measured alternatives after the first working vertical slice.
- Keep GPU targets optional so the portable CPU build never requires a device SDK.
- Upload immutable scene buffers once, keep path state, accumulation, and any later work queues
  resident, and transfer only control data and completed film layers during a render.
- Start with triangle traversal, a minimal opaque material, and environment lighting; add textures,
  emissive triangles, transmission, cutouts, and post-processing only after parity at each step.

Status: implemented for the complete current CPU estimator on the Windows AMD development device.
The focused G0 through G3b probes remain available, while `VulkanPathRenderer` adds tiled SPP
accumulation, complete material textures, candidate alpha cutouts, environment and area-light
sampling, emission, reflection, transmission, absorption, nested media, and a 16-bounce iterative
path. It returns the common Film representation and the `raym0nade_gpu_render` application feeds
that result into the existing host post-processing path. Cross-platform GPU validation,
stored-image acceptance thresholds, proactive memory-budget policy, and discrete-GPU performance
remain open parts of the stage gate.

Gate: every supported feature has a CPU/GPU comparison scene, and unsupported features fail through
capability checks rather than silently changing appearance.

### Stage 5: Optimize and specialize

- Evaluate SAH/wide BVHs, hardware ray-tracing APIs, texture compression/samplers, persistent
  kernels, queue compaction, adaptive sampling, and device-side denoising from profiles.
- Add backend-specific caches keyed by packed scene version without changing application ownership.
- Track image error, rays/s, build/upload time, memory, and energy use; do not accept performance
  changes based only on wall-clock anecdotes.

Gate: optimized backends remain within the documented reference tolerance and retain the CPU
fallback on every supported host platform.

## Rules for extending the architecture

- Keep repository source and documentation in English.
- Add public APIs only under `include/raym0nade/`; keep Assimp, Python, libpng, platform, and future
  GPU SDK headers out of portable public headers whenever possible.
- Express ownership with values and RAII. When borrowing is necessary, document the lifetime and
  prefer indices/spans over raw pointers for new scene structures.
- Validate at subsystem boundaries, then keep hot inner loops simple and measurable.
- Preserve deterministic random-stream keys when changing scheduling or parallelism.
- Keep CPU-only configuration buildable and tested when optional backend dependencies are absent.
- Update this document when a module boundary or lifetime contract changes, and update
  `development-log.md` when a phase, validation result, or known limitation changes.
