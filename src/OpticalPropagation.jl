module OpticalPropagation

using FFTW

export
    angular_spectrum,
    fresnel_diffraction_double,
    fresnel_diffraction_single

include("angular_spectrum.jl")
include("fresnel_diffraction.jl")

end # module
