# Repository Agent Guidelines

These instructions apply to the entire repository and are intended to survive long-running
maintenance and context compaction.

## Communication and repository language

- Communicate with the project owner in Chinese unless they request another language.
- Keep all repository content written or rewritten by agents in English, including source,
  comments, documentation, scripts, test names, commit messages, and generated configuration.
- Record chronological progress, validation results, benchmarks, decisions, and open risks in
  `docs/development-log.md`. Keep durable architecture contracts in `docs/architecture.md`.

## Branch and change discipline

- Modernization work belongs on the `dev` branch unless the project owner selects another branch.
- Preserve unrelated user changes and large local assets. The `model/`, `output/`, and `build/`
  directories are working data and must not be committed.
- Keep `scripts/` limited to runtime decoder modules required by the renderer. Use Conda and CMake
  commands directly instead of adding workstation-specific setup, build, or launcher wrappers.
- Do not commit generated binaries, rendered images, downloaded archives, dependency caches, or
  machine-specific IDE files.
- Prefer small, reviewable changes with tests. Do not silently preserve known undefined behavior
  merely to match legacy implementation details.

## Portability and dependencies

- Keep supported source and CMake behavior portable across Windows, Linux, and macOS.
- Use standard C++17 and avoid platform-specific APIs unless they are isolated behind a guarded
  implementation with a portable alternative.
- Manage Assimp, libpng, Python, and image-decoder packages through the `raym0nade` Conda
  environment while they remain dependencies. Keep third-party types out of public headers.
- Dependency downloads must use a direct connection. Setup/build processes must not route package
  traffic through the local `127.0.0.1:7897` proxy and must not make persistent proxy changes.
- CMake must not download dependencies during configuration. Keep build products under `build/`.

## Renderer correctness and performance

- Treat the CPU renderer as the correctness reference for future GPU backends.
- Sampling must never rely on random endpoints, unchecked zero denominators, non-finite values, or
  invalid probability distributions. Estimator changes require focused tests and visual regression.
- A fixed scene, settings object, and seed must produce pixel-identical output regardless of CPU
  worker count. Preserve row/coordinate-based seeding when changing scheduling.
- Keep `Model` immutable after construction. Its faces borrow stable material and vertex-data
  addresses, and its BVH borrows the final face array; never reallocate those arrays after BVH
  construction.
- Measure performance changes with the same scene, camera, resolution, sample count, seed, and
  build type. Report renderer time separately from scene loading and post-processing.
- Do not begin a GPU backend by duplicating the current pointer-rich CPU scene. First preserve a
  tested CPU baseline and introduce indexed, packed, backend-neutral scene data.

## Required validation

- Configure, build, and run CTest in both Debug and Release for changes that affect C++ or CMake.
- Keep warning levels enabled (`/W4 /permissive-` on MSVC and `-Wall -Wextra -Wpedantic` on GCC or
  Clang) and resolve new warnings.
- Run `git diff --check` and scan authored files for accidental non-English text before handoff.
- For rendering or sampling changes, run a small deterministic 1-thread versus multi-thread image
  comparison and a representative Bistro smoke render when the local asset is available.
- Update `docs/development-log.md` only after the corresponding validation has actually completed.

## Common commands

```text
conda activate raym0nade
cmake --workflow --preset debug
cmake --workflow --preset release
```

On Windows, run the workflow commands from an x64 Native Tools Command Prompt for Visual Studio
2022. Run executables while the Conda environment remains active.

Read `README.md`, `docs/architecture.md`, and `docs/development-log.md` before continuing a
partially completed modernization phase.
