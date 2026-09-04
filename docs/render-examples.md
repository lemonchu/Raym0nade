# Render Recipes and Legacy Records

This document separates runnable settings for the current renderer from values preserved only as
historical records. Current recipes are intended to reproduce the overall appearance of the old
project at a useful modern preview size. They are not bitwise or radiometrically equivalent to the
legacy estimator.

## Interactive input contract

The current command sequence is:

```text
create model <model-id>
create settings <settings-id>
render <model-id> <settings-id>
exit
```

`create model` reads the scene directory, scene filename, and HDR filename in that order. Use
`null` instead of an HDR filename for a scene without an environment.

`create settings` then reads:

1. Camera direction, right, and up vectors.
2. Camera-position coefficients along those three vectors.
3. Pixel scale, focus distance, circle-of-confusion radius, and exposure.
4. Width and height.
5. Samples per pixel, worker threads, and direct-light probability.
6. Output path prefix.

Paths and identifiers are whitespace-delimited. Use `0` worker threads to select the available
hardware concurrency. Exposure is a linear multiplier before highlight compression and display
gamma; it is not an EV value. Samples per pixel and direct-light probability control estimator
allocation and variance, not intended brightness.

## Resolution and depth-of-field scaling

Pixel scale is an image-plane step per pixel, and circle of confusion is a radius measured in
pixels. For an aspect-ratio-preserving resize, let `s = newWidth / oldWidth`:

```text
newPixelScale = oldPixelScale / s
newCircleOfConfusion = oldCircleOfConfusion * s
```

Focus distance and exposure do not change merely because resolution changes. For example, the
2048 x 1152 exterior camera uses pixel scale `0.001025` and CoC `32`; the corresponding 1024 x 576
values are `0.00205` and `16`, while the 512 x 288 values are `0.0041` and `8`.

The daylight recipes below deliberately set CoC to zero while exposure is evaluated. Depth of
field and bloom spread bright pixels and make exposure comparisons harder.

## Choosing an output image

Use `<prefix>(Filter_FXAA).png` to evaluate exposure and color. Bloom variants deliberately spread
highlights and should not be used to judge base exposure. `Raw` still includes the renderer's
spatial outlier clamp and is not an untouched estimator dump. Depth-of-field variants are emitted
only when CoC is positive.

Before comparing Bistro images, confirm that the loader reports the complete known topology:
8,496,360 vertices and 2,832,120 faces. The installed Assimp 5.4.3 can intermittently omit complete
meshes with both the legacy and current importer settings. A run with different counts is not a
valid appearance or timing baseline.

## Current daylight appearance recipe

The primary current appearance trial uses exposure `12`. Exposure `20` was rejected as visibly too
bright. Exposure `14` is an owner-approved, slightly brighter alternative. The move to `12` is the
next requested darker trial, with all other camera and sampling values held constant.

### 512 x 288 preview

The complete input is also stored in `examples/bistro_daylight_appearance_512.txt`.

```text
create model bistro_daylight
model/Bistro_v5_2/
BistroExterior.fbx
san_giuseppe_bridge_4k.hdr
create settings daylight_appearance_512
0.96718 -0.2 -0.156768
0.16 0 0.987117
-0.1974234 -0.9798 0.032
-16.5 -0.5 -2.25
0.0041 15.0 0.0 12.0
512 288
64 0 0.7
output/BistroDaylightAppearance512
render bistro_daylight daylight_appearance_512
exit
```

The validated preview imported 8,496,360 vertices and 2,832,120 faces. It used 7,896,858 direct
light samples, took 14.546 seconds for the render kernel and 16.409 seconds including export on the
calibration machine, and produced the following `Filter_FXAA` display-luminance statistics:

| Metric | Exposure 12 preview |
| --- | ---: |
| Mean | 0.507618 |
| Median | 0.490196 |
| 75th percentile | 0.656863 |
| 90th percentile | 0.927843 |
| Pixels with at least one clipped channel | 10.272% |

### 1024 x 576 quality recipe

The complete input is also stored in `examples/bistro_daylight_appearance_1024.txt`.

```text
create model bistro_daylight
model/Bistro_v5_2/
BistroExterior.fbx
san_giuseppe_bridge_4k.hdr
create settings daylight_appearance_1024
0.96718 -0.2 -0.156768
0.16 0 0.987117
-0.1974234 -0.9798 0.032
-16.5 -0.5 -2.25
0.00205 15.0 0.0 12.0
1024 576
64 0 0.7
output/BistroDaylightAppearance1024
render bistro_daylight daylight_appearance_1024
exit
```

Use exposure `14` in either recipe when a slightly brighter daylight presentation is desired. Do
not change pixel scale when changing exposure.

## Daylight exposure derivation

The old sky sampler stored each HDR texel with the approximate solid angle

```text
sin(phi) * 2*pi / (width*height)
```

while the current renderer uses the exact latitude-longitude texel solid angle. For this 4096 x
2048 HDR, the constant ratio between the exact and legacy midpoint expressions is:

```text
2 * 2048 * sin(pi / (2 * 2048)) = 3.141592346...
50 / 3.141592346...              = 15.915496...
```

Therefore `15.9155` is the precise correction for that solid-angle term alone. It is not the final
appearance exposure. The legacy sampler also retried directions until they entered the shaded
hemisphere without compensating the conditional probability, while the current renderer allows
the HDR environment and emissive mesh lights to contribute together. Both effects vary by surface
normal and image region. They cannot be folded into one exact global multiplier.

The exposure `12` trial is consequently a visual compatibility value: it starts from the exact
solid-angle calculation, leaves headroom for the additional current light contribution, and follows
the requested darker direction after exposure `14`. It must not be described as `50 / pi` or as a
general conversion for other HDR images and scenes.

## Scene definitions

Use different model identifiers when more than one scene is loaded in the same process:

| ID | Scene | Environment |
| --- | --- | --- |
| `bistro_daylight` | `BistroExterior.fbx` | `san_giuseppe_bridge_4k.hdr` |
| `bistro_night` | `BistroExterior.fbx` | `null` |
| `bistro_interior` | `BistroInterior_Wine.fbx` | `null` |

## Archived legacy records

The following values are transcribed from commit `32a4883`. They document the old project and are
not automatically valid current recipes. In particular, the historical exterior overview and
alternate table records contain zero samples per pixel. The old implementation accepted zero but
produced no sampled surface lighting; the current renderer correctly rejects it.

| Legacy ID | Historical scene pairing | Resolution | Pixel scale | Focus | CoC | Exposure | SPP | Direct probability |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `arg_exterior` | Exterior night, `null` | 2048 x 1152 | 0.001025 | 15 | 32 | 200 | 0 | 0.7 |
| `arg_exterior_table` | Exterior night, `null` | 2048 x 1152 | 0.001 | 1.5 | 10 | 36 | 256 | 0.7 |
| `arg_interior` | Interior, `null` | 4096 x 2304 | 0.00032 | 0 | 0 | 512 | 320 | 0.75 |
| `arg_interior_table` | Interior, `null` | 2048 x 1152 | 0.00048 | 0 | 0 | 512 | 320 | 0.7 |
| `arg_interior_counter` | Interior, `null` | 2048 x 1152 | 0.0005 | 0 | 0 | 64 | 640 | 0.7 |
| `arg_table5` | Exterior; environment not recorded | 1024 x 576 | 0.0012 | 16 | 8 | 50 | 0 | 0.7 |

The exposure `200` record is explicitly labelled as a night exterior. It must not be paired with
the daylight HDR and presented as a historical daylight setting. The history contains no complete,
runnable daylight camera/exposure preset; the current daylight recipe above is derived and visually
selected.

## Compatibility limitations

No global parameter set can make every current pixel identical to the legacy output. The current
renderer intentionally fixes invalid PDFs, endpoint-sensitive random sampling, transmission and
reflection sampling, material-factor import, finite-value handling, and environment solid angles.
It also supports environment and mesh-emitter next-event estimation at the same time. These changes
alter local light ratios rather than applying one constant scale.

For future appearance calibration:

1. Keep the scene, camera, resolution, pixel scale, CoC, SPP, and direct probability fixed.
2. Compare `Filter_FXAA`, with bloom and depth of field disabled.
3. Search exposure only after confirming identical imported topology.
4. Record owner-approved settings as scene-specific appearance recipes, not universal legacy
   conversion factors.
