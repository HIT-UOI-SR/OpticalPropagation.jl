"""
    angularspectrum(Uin, d, λ, lx, ly)

calculate the propagation light field based on the angular spectrum.

The input and output use the same uniformly sampled periodic window. Spatial
frequencies follow FFT order for both odd and even sample counts. Nonpropagating
components (`(λ*fx)^2 + (λ*fy)^2 > 1`) are suppressed for nonzero distances,
including backward propagation. Zero distance returns an independent complex copy.
`d` must be finite; `λ`, `lx` and `ly` must be finite and positive.

## Arguments

- `Uin::AbstractArray{<:Number,2}`: Complex array of input complex amplitude.
- `d::Real`: Distance to propagate in metres.
- `λ::Real`: Wavelength of light to propagate.
- `lx::Real`: The physical size of the input data along the x-axis.
- `ly::Real`: The physical size of the input data along the y-axis.

## Returns

- `::Array{<:Number,2}`: Complex amplitude data after propagation.
"""
function angularspectrum(Uin::AbstractArray{<:Number,2}, d::Real, λ::Real, lx::Real, ly::Real)
    n, m = _propagation_grid(Uin, d, λ, lx, ly)
    spectrum = Matrix{ComplexF64}(Uin)
    d == 0 && return spectrum
    fft!(spectrum)
    fx, fy = fftfreq(n, n/lx), fftfreq(m, m/ly)
    carrier = cispi(2*d/λ)
    for j in 1:m, i in 1:n
        transverse = (λ*fx[i])^2 + (λ*fy[j])^2
        if transverse <= 1
            # sqrt(1-s)-1 = -s/(1+sqrt(1-s)) avoids cancellation near DC.
            relative_phase = -2*d/λ * transverse/(1 + sqrt(1 - transverse))
            spectrum[i, j] *= carrier * cispi(relative_phase)
        else
            spectrum[i, j] = 0
        end
    end
    ifft!(spectrum)
end


"""
    angularspectrum(Uin::MonoLightField2D, d::Unitful.Length) -> MonoLightField2D

calculate the propagation light field for [`MonoLightField2D`](@ref) based on the angular spectrum.
"""
function angularspectrum(Uin::MonoLightField2D, d::Unitful.Length)
    Uout = angularspectrum(Uin.data, auval(d), auval(Uin.wavelength), auval(Uin.size[1]), auval(Uin.size[2]))
    _withdata(Uin, Uout)
end

"""
    angularspectrum(d)

represent an [`angularspectrum`](@ref) propagation functor with the specified distance.
"""
angularspectrum(d) = x->angularspectrum(x, d)
