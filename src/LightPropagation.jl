module LightPropagation

using OpticalTransform

export
    angular_spectrum,
    fresnel_diffraction_double

include("angular_spectrum.jl")
include("fresnel_diffraction.jl")

end # module
