"""
A package of useful optical propagation methods and simple 2-dimensional light field data types.
"""
module OpticalPropagation

using FFTW
using Unitful
using RecipesBase

import Base:
    +, -, *, /, \

export
    AbstractLightField2D,
    MonoLightField2D,
    angularspectrum,
    fresnel2,
    fresnel1

include("lightfield2d.jl")
include("angularspectrum.jl")
include("fresnel.jl")
include("utils.jl")

end # module
