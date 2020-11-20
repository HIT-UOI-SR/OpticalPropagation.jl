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
