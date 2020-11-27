# 2D Light Field

The 2-dimensional light field data types provide a more convenient way to represent the light field in the plane space and can be used with our optical propagation methods.

```@docs
    MonoLightField2D
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
