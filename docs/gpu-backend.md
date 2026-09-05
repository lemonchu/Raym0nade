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

The implemented packed format version 2 contains the vertex, triangle-index, triangle-material, and
fixed material arrays. Material flags explicitly preserve whether the source uses alpha cutout or
a diffuse, specular, emissive, or normal texture. All texture IDs deliberately remain invalid
because GPU texture storage and sampling are not implemented yet. Separating presence from IDs
lets capability checks reject unsupported shading instead of silently substituting constants.
Referenced alpha cutouts are rejected at the Vulkan geometry boundary rather than being treated as
opaque.

The private `VulkanRuntime` is common infrastructure for the G2 intersector and G3a renderer. It
owns the selected device, queue synchronization, persistent device-local vertex, index,
triangle-material, and material buffers, and persistent BLAS/TLAS. The scratch and instance-input
buffers used only while building the acceleration structures are released after the build
submission completes.

The CPU BVH will eventually reorder a primitive-ID array instead of geometry. The Vulkan backend
builds its own hardware acceleration structures from the same stable vertex and index arrays; the
current median-split CPU BVH is not uploaded.

Texture deduplication is a first-order requirement. The local Bistro asset contains hundreds of DDS
files and about 1.36 GiB of compressed texture input. GPU upload must query the memory budget,
account for decoded mip storage plus AS scratch space, and fail explicitly rather than oversubscribe
memory silently.

## Integrator migration

Mechanical function-by-function GPU annotation is not viable. The hot path is converted in this
order:

1. Separate integration into `renderToFilm` and keep output naming, post-processing, and PNG I/O on
   the host.
2. Introduce indexed `SceneData` and make the CPU renderer consume it without changing images.
3. Replace dynamic path state with fixed-capacity records. The maximum medium depth follows the
   path-depth limit, and temporary light contributions use a checked fixed capacity.
4. Replace recursive transport with an iterative state machine on the CPU and revalidate it.
5. Replace sequential `std::mt19937` consumption with a counter-based stream keyed by seed, pixel,
   sample, bounce, and dimension.
6. Port ray generation, intersection, minimal shading, and accumulation into one Vulkan compute
   kernel.
7. Add environment and area next-event estimation, textures and normal maps, reflection,
   transmission and nested media, then alpha cutouts in separately validated steps.
8. Profile megakernel execution before deciding whether to split generate, intersect, shade,
   shadow, compact, and accumulate queues.

For alpha cutouts, the Ray Query candidate loop can sample alpha using the candidate primitive and
barycentrics, then confirm or reject the candidate without restarting traversal. Feature parity
must retain the current double-sided geometry, cutoff, maximum skipped-layer policy, texture-level
behavior, and shadow semantics until an intentional estimator change is separately reviewed.

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

GPU `BaseColor` is deliberately restricted to referenced opaque materials without diffuse
textures. `ShapeNormal` supports every scene accepted by the current Vulkan geometry boundary.
Lighting, random sampling, path continuation, texture sampling, accumulation, and post-processing
remain outside G3a.

### G3: Minimal complete render - in progress

- Render opaque diffuse surfaces with a constant or environment light entirely on the GPU.
- Return linear Film/AOV buffers rather than only tone-mapped PNG data.
- Compare against a high-SPP CPU reference with an error tolerance based on CPU seed variance.
- Keep the result repeatable for a fixed GPU, driver, shader binary, settings, and seed.

G3a closes the camera, primary traversal, constant base-color, geometric-normal, and linear-image
boundary only. G3 remains open until lighting and path integration execute on the GPU under a
corresponding CPU correctness comparison.

### G4: Feature parity

- Add textures, mip filtering, emissive triangles, normal maps, reflection, cutouts, transmission,
  absorption, and medium nesting one feature at a time.
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
numbers are host wall-clock diagnostics for tiny ray batches; G3a adds valid GPU dispatch
timestamps, but only for a 4 x 4 primary-AOV diagnostic. None satisfies the performance protocol
above or demonstrates renderer speedup. A future supported AMD discrete GPU uses provisional
targets of at least 5x path-kernel speedup and 3x end-to-end speedup over the 16-core CPU reference.
Missing those targets triggers profiling and design review before more features are ported.
