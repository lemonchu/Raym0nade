# AMD GPU Backend

This document defines the experimental GPU implementation on the `dev-gpu` branch. It records the
selected API, host/device boundary, data contracts, migration decisions, validation gates, and
performance methodology. The CPU renderer remains the correctness reference and the portable
fallback.

## Decision

The first GPU backend is headless Vulkan compute with `VK_KHR_ray_query`, initially restricted to
AMD devices by vendor ID. It requires no window, surface, swapchain, or graphics pipeline. The
implementation uses only KHR functionality so a later policy change can permit other vendors
without rewriting the backend.

Ray Query was selected ahead of a full Ray Tracing Pipeline because it uses the same BLAS and TLAS
while avoiding shader binding tables and multiple shader stages. It also maps naturally to an
iterative per-pixel kernel and a later wavefront queue design. A Ray Tracing Pipeline remains a
measured A/B candidate after the first opaque-material vertical slice; it must not be assumed
faster or slower without profiling.

HIPRT can remove acceleration-structure construction and traversal boilerplate and offers a more
C++-like kernel language, but it does not remove the need to rewrite the pointer-rich scene,
recursive integrator, random generator, textures, lights, medium state, and film accumulation. The
current `gfx1036` APU is not in AMD's official Windows HIP SDK support matrix. HIPRT is therefore
deferred until an officially supported AMD discrete GPU is available for a controlled comparison.

Relevant primary references:

- [Khronos Vulkan ray tracing guide](https://docs.vulkan.org/guide/latest/extensions/ray_tracing.html)
- [Khronos Ray Query reference](https://docs.vulkan.org/refpages/latest/refpages/source/VK_KHR_ray_query.html)
- [AMD HIPRT](https://gpuopen.com/hiprt/)
- [AMD HIP SDK Windows requirements](https://rocm.docs.amd.com/projects/install-on-windows/en/latest/reference/system-requirements.html)
- [AMD Ryzen 9 9950X specifications](https://www.amd.com/en/products/processors/desktops/ryzen/9000-series/amd-ryzen-9-9950x.html)

## Non-negotiable execution boundary

The backend boundary is a complete render, not one ray or one bounce:

```text
Assimp and image decoding on the CPU
                |
                v
   indexed immutable PackedSceneData
                |
        one upload / AS build
                |
                v
GPU ray generation -> traversal -> shading -> path continuation -> Film accumulation
                |
       completed Film tile readbacks
                |
                v
      CPU post-processing and PNG export
```

Scene buffers and acceleration structures stay resident for the lifetime of a renderer. Within a
tile, path state and four radiance accumulators stay resident across every SPP batch. There is no
CPU round trip per ray, visibility query, or bounce.

`VulkanPrimaryRenderer` retains the earlier single-dispatch `LinearImage` diagnostic path.
`VulkanPathRenderer` implements the complete-render boundary with one invocation per pixel, an
iterative 16-bounce loop, and configurable SPP batches. Each invocation is the only writer for its
tile pixel, avoiding floating-point atomics. The first batch initializes the G-buffer and
accumulators; later batches resume them, and only the completed tile is copied to the host Film.

One command buffer records every sample batch for a tile, with shader-write to
shader-read/write barriers between resumed accumulations and a final transfer copy. The tile uses
one queue submission and fence wait. The default batch contains eight samples and the public limit
is 64. This removes ordinary per-batch host waits; it does not yet bound an extreme high-SPP tile
submission. A wavefront pipeline is introduced only after profiles show that megakernel divergence
or register pressure justifies its additional queues, storage, and synchronization.

Path rendering accepts a positive compute-queue count, with one as the portable default. The
runtime first prefers a capacity-sufficient compute-only family, then a graphics-and-compute
family, and never silently lowers the request. Selected queues share the immutable scene, BLAS/TLAS,
and pipeline. Each queue owns its command, synchronization, descriptor, output, readback, and host
tile state; an atomic counter assigns independent tiles. Queue count changes scheduling, not Philox
addresses or Film values.

## Backend-neutral scene contract

The existing `Face` graph cannot cross the device boundary: it contains borrowed `VertexData*` and
`Material*` pointers, the CPU BVH reorders faces, and emissive lights copy faces carrying those
pointers. The replacement scene uses indexed containers with trivially-copyable GPU-visible
records. A backend treats a validated packed scene as immutable for its upload lifetime.

The initial layout is:

| Array | Contract |
| --- | --- |
| Vertices | 32-byte, 16-byte-aligned records stored as two float4 values: `(position.xyz, normal.x)` and `(normal.yz, uv.xy)`. The same buffer can be consumed by the Vulkan AS builder with a 32-byte vertex stride. |
| Triangle indices | Tightly packed `uint32_t[3 * triangleCount]`, as required by the AS triangle input. |
| Triangle materials | One `uint32_t` material ID per triangle. |
| Materials | Fixed-size records containing diffuse/opacity, emission/IOR, transmission/roughness, metallic/specular/flags, and four texture IDs. |
| Textures | A shared `TextureStore`, deduplicated by normalized source path. Materials contain IDs rather than owning copies. |
| Area lights | Light records plus global emissive-triangle ID and normalized CDF arrays; no copied triangles. |
| Environment | Linear radiance, normalized importance CDF, dimensions, and row solid angles. |

All IDs use `uint32_t`; `0xffffffff` is the invalid sentinel. Import validates finite values, index
ranges, allocation overflow, and count limits before publishing a scene. GPU-visible ABI records
require compile-time checks for size, alignment, standard layout, and trivial copyability.

The implemented packed format version 4 contains vertex, triangle-index, triangle-material, fixed
material, texture, mip, encoded RGBA8 texel, area-light, area-light-triangle, environment-row, and
environment-texel arrays. Textures are deduplicated by normalized source path without folding case.
Every present material texture has a checked ID, and each texture owns one contiguous, complete mip
range whose levels own exact contiguous texel ranges. The encoded-byte representation intentionally
preserves the CPU sampler's alpha and gamma-decode behavior.

Area lights refer to packed vertex and emissive-material IDs instead of copied host faces. Their
records preserve center, power, geometric triangle area, normalized face probability, and CDF. The
environment representation stores linear HDR radiance, exact per-row texel solid angle, a row
distribution, and a conditional CDF within each row. Both construction and validation derive the
environment distribution from the CPU luminance function multiplied by exact solid angle. The
validator also recomputes area-light geometry and rejects broken ABI metadata, ranges, reserved
fields, distributions, or non-finite values.

The private `VulkanRuntime` is common infrastructure for the arbitrary-ray, primary-AOV, and
path-rendering APIs.
It owns the selected device, queue synchronization, persistent device-local vertex, index,
triangle-material, material, texture, light, and environment buffers, and persistent BLAS/TLAS.
One reusable host-visible staging allocation transfers every array in chunks no larger than 16 MiB,
so decoded texture storage does not require an equally large staging duplicate. The scratch and
instance-input buffers used only while building the acceleration structures are released after the
build submission completes.

The logical encoded-texel array uses buffer-device-address paging rather than one large descriptor.
The runtime chooses a power-of-two page size no larger than the requested size (256 MiB by default)
or `maxStorageBufferRange`, allocates device-local pages, and uploads a compact page table with
64-bit device addresses. Shader lookup uses the stored page shift and mask to resolve a logical
`uint32_t` texel index through `GL_EXT_buffer_reference`. The runtime validates table size,
address alignment, page ranges, total `uint32_t` indexability, and
`maxMemoryAllocationCount`. All other storage buffers remain individually checked against
`maxStorageBufferRange`.

The path and primary shaders reuse texture descriptors, wrapped UVs, mip decisions, and paging
metadata within each lookup. A zero mip blend skips the second bilinear sample. Four same-page
bilinear texels share one resolved page address, while cross-page taps keep the exact per-texel
fallback.

Mixed opaque/cutout scenes build two geometries in one BLAS. Opaque triangles use
`VK_GEOMETRY_OPAQUE_BIT_KHR`; cutout triangles retain candidate confirmation with
`VK_GEOMETRY_NO_DUPLICATE_ANY_HIT_INVOCATION_BIT_KHR`. Geometry-local primitive IDs are remapped
to original packed triangle IDs in the shader. All-opaque and all-cutout scenes use one geometry;
the old unified candidate layout is retained only for explicit benchmark A/B.

The CPU BVH will eventually reorder a primitive-ID array instead of geometry. The Vulkan backend
builds its own hardware acceleration structures from the same stable vertex and index arrays; the
current median-split CPU BVH is not uploaded.

Texture deduplication remains a first-order requirement. Buffer-device-address pages remove the
per-storage-buffer barrier that previously prevented a multi-gigabyte decoded Bistro texture set
from reaching the shader. They do not solve total residency: complete RGBA8 mip chains still exist
in host packed data during setup and in device-local pages afterward, and there is no texture
streaming, eviction, native compressed `VkImage` storage, or proactive heap-budget policy. Native
image sampling is not a mechanical replacement: the CPU mip chain ceil-halves odd dimensions,
whereas Vulkan image mip extents floor-halve them. An exact first experiment therefore needs
single-level integer images per packed mip (optionally array-grouped by equal extent), descriptor
indexing and buffer-device-address fallback variants, suballocated memory, capability gates, and
measured exact-output A/B.

## Integrator migration

Mechanical function-by-function GPU annotation was not viable. The implemented migration is:

1. Separate complete CPU integration into `renderToFilm` and keep output naming,
   post-processing, and PNG I/O in the consuming `exportFilmToFiles` host operation.
2. Derive immutable indexed `PackedSceneData` from `Model`. The CPU renderer remains the
   pointer-rich recursive reference; the Vulkan backend consumes the packed value representation.
3. Replace dynamic device path state with fixed-capacity records. The shader uses a 17-entry medium
   stack for the outside medium plus a maximum of 16 nested path boundaries, fixed random
   dimensions, and four per-pixel radiance/second-moment accumulators.
4. Port the CPU recursion into one iterative GLSL loop through depth 16. This preserves the
   estimator behavior without first rewriting the working CPU reference.
5. Replace sequential `std::mt19937` consumption on the device with Philox4x32-10 addresses keyed
   by seed, pixel, sample, bounce, replicate, and dimension. The CPU and GLSL contracts define the
   same raw integers and endpoint-safe open-(0,1) float conversion independently of scheduling.
6. Execute primary ray generation, Ray Query traversal, shading, continuation, next-event
   visibility, and Film accumulation inside one per-pixel compute megakernel.
7. Consume the complete packed shading boundary: diffuse, specular, emissive, and normal textures;
   candidate alpha cutouts; environment and area-light importance sampling; emission; diffuse,
   microfacet reflection, clearcoat, and transmission response; nested IOR and absorption state;
   direct/indirect mixture compensation; transparent-first-hit replication; roughness
   regularization; and throughput limiting.
8. Read back completed Film tiles and use `finalizeRadianceData` plus the existing CPU
   post-processing pipeline. Profile the megakernel before deciding whether to split generate,
   intersect, shade, shadow, compact, and accumulate queues.

For alpha cutouts, the implemented Ray Query candidate loop samples base-mip alpha using the
candidate primitive and barycentrics, then confirms or rejects the candidate without restarting
traversal. It uses the CPU cutoff of `1e-4` and applies to primary and shadow rays. Vulkan does not
guarantee candidate order, so counting candidates and stopping at 32 could turn an ordinary ray
into a false miss before a nearer opaque hit is visited. The GPU therefore rejects all transparent
candidates and has no layer cap. This is more robust for ordinary traversal but deliberately
differs from the CPU only for a pathological ray that crosses more than 32 distance-ordered
transparent layers: the CPU returns a miss after its defensive cap, while the GPU continues.

## Milestones and gates

### G0: Toolchain and device capability

- CPU presets configure without requiring Vulkan.
- GPU presets find the Conda-managed Vulkan loader, headers, and shader compiler.
- A project-owned probe reports AMD identity, a compute queue, buffer device address,
  acceleration-structure support, and Ray Query support.

Status: complete locally in clean Windows Debug and Release GPU builds on the `0x1002/0x13c0`
integrated GPU. The optional Vulkan targets remain disabled in CPU presets. Cross-platform GPU
compile and runtime gates remain open.

### G1: Hardware intersection

- Compile a checked-in compute shader to SPIR-V under `build/`; do not commit generated SPIR-V.
- Build one BLAS and one identity TLAS.
- Trace deterministic hit and miss rays through one triangle.
- Compare primitive ID, distance, and barycentrics with CPU expectations.
- Run Debug with the Khronos validation layer, explicitly enable synchronization validation, and
  report zero validation errors.

Status: completed locally in clean Windows Debug and Release builds on the `0x1002/0x13c0` AMD
device. Both manual runs enabled validation and synchronization validation and reported zero
errors and zero warnings. Both checked-in compute shaders compile for Vulkan 1.2 and pass
`spirv-val`. Cross-platform GPU compile and runtime gates remain open.

### G2: Packed imported geometry

- Convert the checked-in tiny scene to indexed SceneData.
- Compare CPU BVH and Vulkan primary hits over a deterministic ray set.
- Measure pack, upload, BLAS/TLAS build, dispatch, and readback independently.
- Reject invalid indices, non-finite geometry, unsupported counts, and inadequate memory budgets.

Status: the initial imported-geometry slice is complete locally in clean Windows Debug and Release
builds. A three-material, six-vertex, two-face packed OBJ fixture matched CPU BVH primitive IDs,
distances, and barycentrics across two repeated five-ray batches; validation and synchronization
validation were enabled and reported zero errors and warnings:

| Build | Pack | Host setup | Upload | BLAS/TLAS build | Batch 1 | Batch 2 | Exit |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Debug | 0.019 ms | 81.312 ms | 2.143 ms | 0.882 ms | 0.201 ms | 0.104 ms | 0 |
| Release | 0.003 ms | 78.993 ms | 2.285 ms | 1.040 ms | 0.197 ms | 0.135 ms | 0 |

After this timing-only test update, all five CTest entries passed again in both `gpu-debug` and
`gpu-release`.

This is a functionality slice, not the complete G2 scale gate. A fixture with more than ten faces
must exercise CPU BVH reordering, proactive memory-budget checks remain unimplemented, and exact
ray-interval endpoints and teardown-time validation still need coverage. A fence timeout can also
wait indefinitely while the implementation safely quiesces submitted work.

### G3a: Deterministic primary-AOV render

- Share backend-neutral `ImageExtent`, `PinholeCamera`, `PrimaryAov`, `PrimaryRenderRequest`, and
  row-major `LinearImage` types between CPU and Vulkan.
- Preserve the legacy camera formula exactly: integer pixel locations, no half-pixel offset, and no
  implicit normalization of the supplied camera basis.
- Keep a no-file, no-random CPU `BaseColor` and `ShapeNormal` implementation as the oracle.
- Generate primary rays on the device in one two-dimensional compute dispatch using 8 x 8
  workgroups and read back the completed linear image once.
- Reuse persistent scene buffers and BLAS/TLAS through the private `VulkanRuntime`.
- Report host dispatch/readback time and, when supported, GPU dispatch timestamps with wrap handling
  based on the queue family's timestamp valid-bit count.

Status: complete locally in Windows Debug and Release on `AMD Radeon(TM) Graphics`. CPU Debug and
Release each passed all four CTest entries. GPU Debug and Release each passed all seven entries,
including CPU/GPU `BaseColor` and `ShapeNormal` comparison and repeated-render determinism. The
Release validation run enabled the Khronos layer and synchronization validation and reported zero
errors and zero warnings.

The selected Release diagnostics for the 4 x 4 fixture were:

| Operation | Host dispatch/readback | GPU dispatch timestamp |
| --- | ---: | ---: |
| `BaseColor`, first render | 0.270 ms | 0.009 ms |
| `BaseColor`, second render | 0.104 ms | 0.005 ms |
| `ShapeNormal`, first render | 0.093 ms | 0.006 ms |
| `ShapeNormal`, second render | 0.108 ms | 0.006 ms |

The persistent renderer setup reported 2.861 ms for upload and 0.853 ms for BLAS/TLAS build. These
numbers come from a tiny 4 x 4 correctness diagnostic on the two-compute-unit integrated GPU. They
are not throughput measurements and are not evidence of GPU speedup.

At the G3a checkpoint, GPU `BaseColor` was deliberately restricted to referenced opaque materials
without diffuse textures, while `ShapeNormal` supported the wider geometry boundary. The later
device-ready foundation work below supersedes that temporary restriction with diffuse sampling and
candidate alpha testing. Lighting, random sampling, path continuation, accumulation, and
post-processing remained outside G3a.

### G3b: Deterministic directional DirectDiffuse

- Add `DirectionalLight(directionToLight, incidentRadiance)` to the backend-neutral request and
  add `PrimaryAov::DirectDiffuse` to the CPU and Vulkan implementations.
- Normalize `directionToLight` on the host and use a camera-facing shape normal at the first opaque
  hit. Misses return black.
- Evaluate the deterministic Lambert diagnostic:

  ```text
  max(baseColor, 0) * incidentRadiance * max(dot(N, L), 0) / pi
  ```

- Reuse one `rayQueryEXT` variable inside each invocation: complete the primary query, then issue a
  conditional opaque terminate-on-first-hit shadow query from the hit point.
- Preserve one 8 x 8 two-dimensional dispatch, one submission, and one row-major image readback.
- Carry the request-local light in the 112-byte push-constant block without changing the packed
  scene ABI or descriptor bindings.
- Compute `std::nextafter(kRayEpsilon, +infinity)` on the host and pass it through the padded
  `cameraDirection.w` member as both primary and shadow `tMin`. Vulkan's closed lower bound then
  matches the CPU intersector's strict `t > kRayEpsilon` acceptance rule without changing the
  112-byte push-constant ABI.
- At the G3b checkpoint, reject referenced diffuse textures and non-opaque materials for
  `DirectDiffuse` and reject alpha cutouts rather than treating them as opaque. The subsequent
  device-ready foundation slice below replaces this temporary cutout/texture gate.

Status: complete locally in Windows Debug and Release on `AMD Radeon(TM) Graphics`, based on
checkpoint commit `448d377` (`Add Vulkan GPU rendering foundation`). CPU Debug and Release each
passed all four CTest entries; GPU Debug and Release each passed all seven. The Release direct run
enabled the Khronos validation layer and synchronization validation and reported zero errors and
zero warnings.

After the endpoint-parity correction from final cross-review, the targeted directional-light CTest
selection passed its single test in both GPU Debug and GPU Release (1/1 each). A validation and
synchronization-validation run again reported zero errors and zero warnings. This closes the fixed
G3b primary/shadow lower-bound mismatch only; G2 arbitrary `tMax` and general interval-endpoint
coverage remain open.

The new directional-light fixture and CPU/GPU tests cover analytic hit values, hard shadow,
backlighting, two-sided viewing, direction-scale invariance, a 13 x 9 extent that is not a multiple
of the 8 x 8 workgroup, render-target growth and shrinkage, repeated-render identity, and capability
rejection. The selected 13 x 9 Release diagnostics were:

| Invocation | Host dispatch/readback | GPU dispatch timestamp |
| --- | ---: | ---: |
| `DirectDiffuse`, first | 0.1791 ms | 0.00788 ms |
| `DirectDiffuse`, repeated | 0.1068 ms | 0.00728 ms |

Other 4 x 4 and 13 x 9 diagnostic cases ranged from 0.0903 to 0.2214 ms on the host and from
0.00552 to 0.01188 ms in the GPU timestamp. These tiny correctness workloads do not satisfy the
measurement policy and are not evidence of speedup.

At its checkpoint, G3b deliberately excluded environment and area lighting, emission, specular
response, metallic and roughness response, smooth normals, normal maps, distance attenuation,
random sampling, path continuation, textures, accumulation, and post-processing. It was the first
deterministic GPU lighting-and-shadow slice, not the complete G3 renderer.

### Complete Vulkan path renderer

The foundations after G3b are now consumed by `VulkanPathRenderer`:

- Packed format version 4 supplies deduplicated encoded RGBA8 textures with complete mip chains,
  indexed area lights and triangles, and linear HDR radiance with a hierarchical CDF.
- `VulkanRuntime` uploads texture metadata, buffer-device-address texel pages, lights, and
  environment arrays through a reusable staging buffer capped at 16 MiB per transfer.
- `shaders/include/packed_scene.glsl` centralizes the packed GPU ABI, paged texture sampler, and
  candidate-confirm traversal for both diagnostic and path shaders.
- Material evaluation flips V, repeat-wraps both axes, applies the CPU mip/filter and gamma
  conventions, evaluates all four texture slots, and perturbs the interpolated normal when a
  normal map is present. Cutout traversal uses base-mip alpha and the `1e-4` threshold.
- The path shader consumes area-light and HDR-environment importance distributions for next-event
  estimation and evaluates the same current BSDF, reflection/transmission, medium, regularization,
  and accumulation rules as the CPU renderer.
- Philox addresses are invariant under tile and batch partition changes. Each output pixel stores a
  first-hit G-buffer plus direct/indirect diffuse/specular RGB sums and their second moments.
- Philox blocks are generated once and consumed from a cached four-word block where adjacent
  dimensions share the same counter. Texture lookup reuses descriptor, mip, UV, and page
  calculations, and opaque/cutout triangles use separate BLAS geometries.
- All batches and the readback copy for one ordinary tile use one command buffer submission.
  Optional additional compute queues take independent tiles without changing random addresses.
- Host readback assembles a common `Film`; `finalizeRadianceData` and the existing export path
  provide the same exposure, variance, filtering, display, bloom, depth-of-field, FXAA, and PNG
  behavior as CPU output.
- The compact arbitrary-ray `VulkanRayQueryIntersector` still rejects cutout scenes because its
  public records carry no texture or material-sampling contract; this does not restrict either
  renderer.

Extremely small float probabilities can make adjacent packed CDF entries equal. Shader selection
uses adjacent differences and skips zero-width bins, but the lost probability precision cannot be
recovered after packing.

### G3: Minimal complete render - implementation complete

- True iterative path continuation and stochastic lighting execute on the GPU.
- The renderer returns the common linear Film/G-buffer representation.
- CPU comparison uses statistical tolerance because CPU and GPU intentionally consume different
  random streams.
- A fixed GPU, driver, shader binary, settings object, and seed produce repeatable results; tile and
  sample-batch partitioning do not change Philox addresses.

G3a and G3b remain focused deterministic diagnostics. Full G3 is implemented by
`VulkanPathRenderer`, which performs general current-estimator lighting, multi-bounce transport,
SPP accumulation, and Film readback. The `raym0nade_gpu_render` CLI exposes that path without
adding interactive console policy to the backend.

### G4: Current-estimator feature parity - implementation complete

The GPU path evaluates diffuse, specular, emissive, and normal texture inputs; alpha-tested
traversal; area and environment lights; emission; diffuse and specular response; clearcoat;
reflection; transmission; absorption; nested media; direct/indirect allocation; transparent
replication; and the same 16-bounce limit and preview-oriented clamps as the CPU reference. This
means parity with the repository's current estimator, not with a general Disney or glTF renderer.

Post-processing remains a shared CPU stage after either backend returns Film data. The existing CPU
renderer and all CPU-only builds remain available. Broader analytic/golden coverage, Linux and
macOS GPU validation, memory-budget policy, and discrete-GPU performance are still required before
the experimental backend can be treated as production-ready.

### G5: Performance specialization

- Continue comparing the per-pixel megakernel, a measured wavefront implementation, and an iterative Vulkan
  Ray Tracing Pipeline using the same packed scene and ray workload.
- Use Radeon GPU Profiler for occupancy and stalls, Radeon Raytracing Analyzer for AS/traversal,
  Radeon Memory Visualizer for budgets, and Radeon GPU Analyzer for shader ISA/register pressure.
- Move post-processing to GPU only after path tracing no longer dominates end-to-end time.

Initial accepted specialization now includes Philox block reuse, split opaque/cutout BLAS
geometries, exact paged-texture lookup reuse, one submission per tile, and configurable
multi-queue scheduling. A GPU area-light total-weight cache and a CPU binned-SAH experiment were
reverted after negative measurements. Native Vulkan image sampling and a wavefront path remain
future measured candidates.

## Measurement policy

Every performance report separates:

- import and decode;
- scene packing;
- host-to-device upload;
- BLAS and TLAS construction;
- warm-up;
- GPU kernel time from timestamp queries;
- readback;
- CPU post-processing; and
- image export.

Release measurements disable validation, warm up the workload, run at least ten measured
iterations where practical, and report median, p95, rays/s, peak memory, and image error. Every
Bistro result must record its imported vertex and face counts. The known Assimp topology variation
also affects the legacy CPU path, so it is tracked as a separate importer defect rather than used
to block GPU functionality work. A CPU/GPU timing or image comparison is valid only when both runs
import the same topology.

The current 9950X integrated GPU has only two compute units and is a functionality gate. The G1/G2
numbers are host wall-clock diagnostics for tiny ray batches; G3a/G3b add valid GPU dispatch
timestamps, but only for small correctness diagnostics. A matched Bistro path-render comparison on
this APU has not shown acceleration over the CPU reference. A future supported AMD discrete GPU
uses provisional targets of at least 5x path-kernel speedup and 3x end-to-end speedup over the
16-core CPU reference. Missing those targets triggers profiling and design review before
specialization work is accepted.

The optional `raym0nade_gpu_primary_benchmark` executable supplies a practical wall-clock
harness for the existing `ShapeNormal` slice. It times CPU one-worker, CPU automatic-worker, and
the complete GPU render call outside their public APIs; GPU dispatch/readback and timestamp values
remain separate diagnostics, and PNG export is outside every timed region. It defaults to ten
measured iterations and reports median and nearest-rank p95 statistics. A validation-enabled run
fails if the requested layer is unavailable or reports any errors or warnings. Current Bistro runs
keep the version-4 cutout flags intact, so CPU and GPU both apply diffuse-alpha candidate semantics.
It remains a traversal/AOV harness rather than a path-render benchmark.

The `raym0nade_gpu_render` CLI is the practical beauty-render entry point. It reports import,
packing, scene upload, acceleration-structure build, host render wall-clock, GPU dispatch
timestamps when available, dispatch count, direct-light events, validation state, and total
wall-clock. It intentionally performs one requested render rather than the repeated warm-up and
distribution measurements required by the formal benchmark policy. Default export writes the
filtered/FXAA beauty image; `--all-passes` uses the complete shared Film exporter.

The optional `raym0nade_gpu_path_benchmark` executable supplies a controlled beauty benchmark
harness, but it does not yet implement every item in the formal policy above. It imports one
`Model`, performs fixed-seed CPU warm-ups and repeated renders, packs exactly that topology,
releases the `Model`, and reuses a persistent `VulkanPathRenderer` within each GPU arm. Its report
separates external render-call wall clock, internal host time, GPU timestamps, import, packing,
upload, acceleration-structure construction, and PNG export; it also compares composite
floating-point linear radiance before display processing. The default one warm-up and three
measurements are smoke settings, while performance runs require at least ten measurements.

The current host interval includes command recording, submission, waits, readback, and Film
assembly. Readback is not isolated, and the harness does not measure peak memory or rays per second.
A comparison shader runs from the same packed scene without keeping two device scenes alive.
Primary and comparison arms have independent unified-geometry switches, allowing a split-layout
primary and legacy unified-layout comparison in one import. Because fixed A/B order can be biased
by driver warm-up, DVFS, and thermal state, performance conclusions require both primary-comparison
and `--comparison-first` runs with matching topology.

`--comparison-gpu-queues N` also creates a same-shader comparison arm; unless explicitly
overridden, it inherits the primary geometry layout. The report records each shader's byte size and
a labeled noncryptographic FNV-1a-64 identity. For multiple queues, device timestamps are summed
queue-busy intervals and are not elapsed GPU wall-clock time; external render-call wall clock is
the comparison metric.

On the development APU's complete-topology Bistro workload (512 x 288, 64 SPP, seed 0, exposure 12,
128 x 128 tiles, 64-SPP batches, two warm-ups, ten measurements), matched AB and BA runs found two
queues slower than one by 1.4% and 2.6%, respectively, with bit-identical Film values. The default
therefore remains one, while the positive queue-count option remains available for measurement on
other devices. Two separate same-import split-versus-unified observations found split geometry
faster by 5.7% at 2,800,258 faces and 13.8% at 2,832,120 faces. Because those process topologies
differed, they are supporting observations rather than one paired AB/BA range. Matched-topology
texture AB/BA runs found the exact sampler common-expression changes faster by 2.3% and 4.1%.
These observations are device-specific, not portable performance promises.

The older recorded Full HD result is superseded for image comparison. That run cleared 13 cutout
material flags only in the benchmark-local GPU scene while the CPU continued to see through alpha
cards. CPU and GPU therefore shaded different surfaces around foliage, planters, and windows,
producing maximum absolute error 1.0, mean absolute error 0.007328, RMSE 0.057768, and 4.100% of
pixels outside tolerance. Its error PNG multiplied linear absolute error by 10,000 for visibility,
so a linear difference of only `1e-4` appeared near full brightness. Those numbers diagnose the
old semantic mismatch; they are not an accuracy result for the current candidate-cutout path.

The `raym0nade_gpu_primary_benchmark --gpu-only` mode keeps model import, packed-scene conversion,
Vulkan setup, GPU warm-up, and GPU measurements, but it executes no CPU cold, warm-up, or measured
render. It writes only
`gpu-shape-normal.png`, a GPU-only `timings.csv`, and a summary that contains no CPU comparison or
ratio. This mode is intended for local manual rendering and profiling automation; it does not
change the primary benchmark's `ShapeNormal` scope. Beauty rendering belongs to
`raym0nade_gpu_render`.
