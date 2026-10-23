# 2D Light Field

The 2-dimensional light field data types provide a more convenient way to represent the light field in the plane space and can be used with our optical propagation methods.

```@docs
    MonoLightField2D
```

## Physical coordinates

The first array dimension is x and the second is y. `size(field)` gives sample
counts; the existing `field.size` property and `physicalsize(field)` give physical
window widths. `gridcenter(field)` is the window centre, which defaults to zero.
It is not the arithmetic mean of the sample positions. For each dimension,

```math
\Delta x = L_x/N_x, \qquad x_i = c_x - L_x/2 + (i-1)\Delta x.
```

This preserves the original half-open window convention, including odd sample
counts and singleton dimensions. Coordinate ranges carry length units and use
constant storage. `axes(field)` continues to provide integer array indices.

```@example coordinates
using OpticalPropagation, Unitful

field = MonoLightField2D(ones(4, 3),
    wavelength=633u"nm", size=(8u"mm", 6u"mm"),
    center=(10u"mm", -3u"mm"))
samplepitch(field)
```

```@example coordinates
x, y = coordinates(field)
collect(x), collect(y)
```

Existing constructors still work without `center`. The copy constructor inherits
the centre unless it is explicitly changed. Replacing `size` changes the window
widths around that centre; replacing `data` with a different matrix shape keeps
the same window and changes the sampling pitches. These operations relabel the
samples; they do not interpolate or resample the amplitude. Public constructors
and `copy(field)` make independent amplitude matrices.

`crop(field, I, J)` instead extracts a rectangular region while preserving its
physical sample positions and pitches. It accepts nonempty contiguous unit ranges
and `:` and copies the selected data. Ordinary array slicing retains its existing
behaviour and returns data without light-field metadata.

```@example coordinates
region = crop(field, 2:3, 2:3)
physicalsize(region), gridcenter(region), coordinates(region)
```

Scalar arithmetic and same-grid broadcasting preserve coordinates. Combining
light fields checks sample counts, wavelength, physical widths and window centres,
including in-place assignments. Centre differences are compared at the sample
pitch scale, so a large absolute centre does not hide a displacement of one pixel.
Unlabelled matrices and vectors can be used as masks with ordinary broadcasting
rules, provided they do not expand the physical grid. No interpolation or
coordinate alignment is performed implicitly.

`parent(field)` returns the actual `Matrix{ComplexF64}` without copying. Use this
for matrix FFTs and plans; such operations return numerical data and do not
construct frequency-coordinate metadata. Mutating the parent also mutates the
light field.

```@docs
    physicalsize
    samplepitch
    gridcenter
    coordinates
    crop
```

## Work with Plots.jl

The 2-dimensional light field data types can work with [Plots.jl](http://docs.juliaplots.org/latest/), so it will be more convenient when we want to visualize a light field. For example,

```@example
using OpticalPropagation # hide
using Unitful
using Plots

gauss = MonoLightField2D(
    [exp(-(x^2+y^2)/6400) for x in -200:199, y in -200:199],
    wavelength=632.8u"nm",
    size=(1u"mm", 1u"mm")
)
plot(gauss)
```
