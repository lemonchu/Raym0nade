# AMD GPU Backend Plan

This document defines the experimental GPU work on the `dev-gpu` branch. It records the selected
API, host/device boundary, data contracts, migration order, validation gates, and performance
methodology. The CPU renderer remains the correctness reference and the portable fallback.

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
      indexed immutable SceneData
                |
        one upload / AS build
                |
                v
GPU ray generation -> traversal -> shading -> path continuation -> Film accumulation
                |
        one completed-film readback
                |
                v
      CPU post-processing and PNG export
```

Scene buffers, acceleration structures, path state, random state, and accumulation stay resident
on the GPU for the render. CPU/GPU synchronization at every ray or bounce is forbidden because it
would make dispatch and transfer overhead dominate useful work.

G3a establishes this boundary for deterministic primary AOVs before path state exists. Its
`VulkanPrimaryRenderer` uploads the immutable scene once, generates every primary ray on the device
in one two-dimensional compute dispatch using 8 x 8 workgroups, and performs one completed-image
readback. It never uploads a host-generated ray array. This is an architectural slice of the
intended complete-render boundary, not yet the complete execution flow shown above.

G3b extends the same single-dispatch boundary with deterministic directional Lambert shading. Each
shader invocation performs its primary query and, for a front-lit hit, reuses the same Ray Query
variable for one opaque terminate-on-first-hit shadow query. It still performs one submission and
one completed-image readback; there is no host round trip between visibility and shading.

The first complete renderer uses one invocation per pixel, an iterative bounce loop, and a small
SPP batch per dispatch. Each invocation is the only writer for its pixel, avoiding floating-point
atomics and making repeatability tractable. Small batches bound kernel duration on Windows. A
wavefront pipeline is introduced only after profiles show that megakernel divergence or register
pressure is a limiting factor.

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

The private `VulkanRuntime` is common infrastructure for the G2 intersector and G3a/G3b renderer.
It owns the selected device, queue synchronization, persistent device-local vertex, index,
triangle-material, material, texture, light, and environment buffers, and persistent BLAS/TLAS.
One reusable host-visible staging allocation transfers every array in chunks no larger than 16 MiB,
so decoded texture storage does not require an equally large staging duplicate. The scratch and
instance-input buffers used only while building the acceleration structures are released after the
build submission completes. Every individual storage buffer is checked against
`maxStorageBufferRange` before allocation.

The CPU BVH will eventually reorder a primitive-ID array instead of geometry. The Vulkan backend
builds its own hardware acceleration structures from the same stable vertex and index arrays; the
current median-split CPU BVH is not uploaded.

Texture deduplication is a first-order requirement. The local Bistro asset contains hundreds of DDS
files and about 1.36 GiB of compressed texture input. GPU upload must query the memory budget,
account for decoded mip storage plus AS scratch space, and fail explicitly rather than oversubscribe
memory silently. The current device probe reports a 31.82 GiB local/shared heap, while a separate
local device-limit query reports a per-storage-buffer limit of 4 GiB minus one byte. A read-only
header audit of 622 local Bistro DDS files computes 7,539,069,640 bytes (7.021 GiB) for their
complete expanded RGBA8 mip chains, without
counting another 11 TGA files. The current single texel SSBO therefore cannot represent a
full-quality Bistro upload on this device: texel paging, buffer device address segmentation, or
native compressed `VkImage` storage is required. This calculation is capacity evidence, not a
successful packing or rendering run. Because a Vulkan BLAS is marked candidate geometry when any
scene material uses a cutout, an all-cutout or mixed large scene also remains a traversal-
performance risk.

## Integrator migration

Mechanical function-by-function GPU annotation is not viable. The hot path is converted in this
order:

1. **Foundation complete:** separate complete CPU integration into `renderToFilm` and keep output
   naming, post-processing, and PNG I/O in the consuming `exportFilmToFiles` host operation.
2. Introduce indexed `SceneData` and make the CPU renderer consume it without changing images.
3. Replace dynamic path state with fixed-capacity records. The maximum medium depth follows the
   path-depth limit, and temporary light contributions use a checked fixed capacity.
4. Replace recursive transport with an iterative state machine on the CPU and revalidate it.
5. Replace sequential `std::mt19937` consumption with a counter-based stream keyed by seed, pixel,
    sample, bounce, and dimension. A shared Philox4x32-10 CPU/GLSL contract and CPU known-answer
    tests are complete. A GPU known-answer dispatch now verifies all 12 renderer addresses as raw
    integer and open-(0,1) float bits and verifies repeated-dispatch identity; integrator migration
    remains open.
6. Port ray generation, intersection, minimal shading, and accumulation into one Vulkan compute
   kernel.
7. Add environment and area next-event estimation, remaining material textures and normal maps,
   reflection, transmission, and nested media in separately validated steps. Packed diffuse
   sampling and candidate alpha-cutout traversal are already available to the primary-AOV shader.
8. Profile megakernel execution before deciding whether to split generate, intersect, shade,
   shadow, compact, and accumulate queues.

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

### Post-G3b device-ready foundations

The current foundation extends G3b without claiming a completed path renderer:

- Packed format version 4 supplies deduplicated encoded RGBA8 textures with complete mip chains,
  indexed area lights and triangles, and linear HDR radiance with a hierarchical CDF.
- `VulkanRuntime` persistently uploads texture descriptors, mips, texels, lights, and environment
  arrays through a reusable staging buffer capped at 16 MiB per transfer.
- `shaders/include/packed_scene.glsl` centralizes the packed GPU ABI, texture sampler, and
  candidate-confirm traversal for both the primary diagnostic and future path shaders.
- Diffuse sampling flips V, repeat-wraps both axes, bilinearly or trilinearly filters encoded
  bytes, then decodes RGB with `pow(value, 2.2)`. Primary-ray differentials select the mip
  footprint. Cutout traversal uses base-mip alpha and the `1e-4` threshold.
- The compact arbitrary-ray `VulkanRayQueryIntersector` still rejects cutout scenes because its
  public records carry no texture or material-sampling contract; this does not restrict
  `VulkanPrimaryRenderer`.
- The Philox device test compares 12 full renderer addresses against CPU raw integer and
  open-(0,1) float bits and confirms repeated-dispatch identity.
- `finalizeRadianceData` is the shared host finishing operation for CPU- or GPU-produced Film
  radiance accumulators, keeping exposure, variance, overflow, and non-finite handling in one
  implementation.

The packed lighting arrays are resident but are not yet consumed by a path integrator. Extremely
small float probabilities can also make adjacent CDF entries equal; shader searches must tolerate
these plateaus.

### G3: Minimal complete render - in progress

- Extend beyond the single directional `DirectDiffuse` diagnostic with constant environment or
  emission, or introduce true iterative path state and continuation entirely on the GPU.
- Return linear Film/AOV buffers rather than only tone-mapped PNG data.
- Compare against a high-SPP CPU reference with an error tolerance based on CPU seed variance.
- Keep the result repeatable for a fixed GPU, driver, shader binary, settings, and seed.

G3a closes the camera, primary traversal, constant base-color, geometric-normal, and linear-image
boundary. G3b closes a deterministic directional Lambert and hard-shadow sub-gate. Full G3 remains
open until broader lighting and path integration execute on the GPU under a corresponding CPU
correctness comparison. Packed textures, candidate cutouts, uploaded lighting data, the no-file
Film boundary, and the verified counter-RNG contract are foundations; none alone constitutes GPU
beauty rendering.

### G4: Feature parity

- Reuse the completed diffuse-texture, mip-filtering, and cutout-traversal foundations while adding
  emissive evaluation, remaining texture slots, normal maps, reflection, transmission, absorption,
  and medium nesting one feature at a time.
- Unsupported features fail capability validation rather than silently changing appearance.
- The existing CPU renderer and every platform's CPU-only build remain usable.

### G5: Performance specialization

- Compare the per-pixel megakernel, a measured wavefront implementation, and an iterative Vulkan
  Ray Tracing Pipeline using the same packed scene and ray workload.
- Use Radeon GPU Profiler for occupancy and stalls, Radeon Raytracing Analyzer for AS/traversal,
  Radeon Memory Visualizer for budgets, and Radeon GPU Analyzer for shader ISA/register pressure.
- Move post-processing to GPU only after path tracing no longer dominates end-to-end time.

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
iterations where practical, and report median, p95, rays/s, peak memory, and image error. Bistro is
accepted only when import reports the complete known topology of 8,496,360 vertices and 2,832,120
faces.

The current 9950X integrated GPU has only two compute units and is a functionality gate. The G1/G2
numbers are host wall-clock diagnostics for tiny ray batches; G3a/G3b add valid GPU dispatch
timestamps, but only for 4 x 4 and 13 x 9 correctness diagnostics. None satisfies the performance
protocol above or demonstrates renderer speedup. A future supported AMD discrete GPU uses
provisional targets of at least 5x path-kernel speedup and 3x end-to-end speedup over the 16-core CPU
reference. Missing those targets triggers profiling and design review before more features are
ported.

The optional `raym0nade_gpu_primary_benchmark` executable supplies a practical wall-clock
harness for the existing `ShapeNormal` slice. It times CPU one-worker, CPU automatic-worker, and
the complete GPU render call outside their public APIs; GPU dispatch/readback and timestamp values
remain separate diagnostics, and PNG export is outside every timed region. It defaults to ten
measured iterations and reports median and nearest-rank p95 statistics. A validation-enabled run
fails if the requested layer is unavailable or reports any errors or warnings. Current Bistro runs
keep the version-4 cutout flags intact, so CPU and GPU both apply diffuse-alpha candidate semantics.
They still do not satisfy the complete-topology or textured-beauty gates above.

The older recorded Full HD result is superseded for image comparison. That run cleared 13 cutout
material flags only in the benchmark-local GPU scene while the CPU continued to see through alpha
cards. CPU and GPU therefore shaded different surfaces around foliage, planters, and windows,
producing maximum absolute error 1.0, mean absolute error 0.007328, RMSE 0.057768, and 4.100% of
pixels outside tolerance. Its error PNG multiplied linear absolute error by 10,000 for visibility,
so a linear difference of only `1e-4` appeared near full brightness. Those numbers diagnose the
old semantic mismatch; they are not an accuracy result for the current candidate-cutout path.

The `--gpu-only` mode keeps model import, packed-scene conversion, Vulkan setup, GPU warm-up, and
GPU measurements, but it executes no CPU cold, warm-up, or measured render. It writes only
`gpu-shape-normal.png`, a GPU-only `timings.csv`, and a summary that contains no CPU comparison or
ratio. This mode is intended for local manual rendering and profiling automation; it does not
expand the supported AOV, material, lighting, or path-integration feature set.
