module OpticalPropagation

using FFTW
using Unitful

import Base:
    +, -, *, /, \

export
    AbstractLightField2D,
    MonoLightField2D,
    angular_spectrum,
    fresnel_diffraction_double,
    fresnel_diffraction_single

include("types.jl")
include("angular_spectrum.jl")
include("fresnel_diffraction.jl")

end # module
