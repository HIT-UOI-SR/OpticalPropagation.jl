# OpticalPropagation.jl

A package of useful optical propagation methods and simple 2-dimensional light field data types.

## Installation

First, add the package

```julia
using Pkg
Pkg.add("OpticalPropagation")
```

To use some light field data types, you need to add [Unitful.jl](https://painterqubits.github.io/Unitful.jl/latest/)

```julia
Pkg.add("Unitful")
```

## Basic Examples

Draw a cross hole and its diffraction pattern:

```@example
using OpticalPropagation
using Unitful
using Plots

cross = MonoLightField2D(
    [Int((abs(x)<10 && abs(y)<50) || (abs(x)<50 && abs(y)<10)) for y in -200:199, x in -200:199],
    wavelength=632.8u"nm",
    size=(1u"mm",1u"mm")
)
patt = angularspectrum(cross, 1u"cm")

plot(
    plot(cross, size=(420,400)),
    plot(patt, size=(420,400)),
    size=(840,400)
)
```
