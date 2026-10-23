"""
    rayleighsommerfeld(Uin, d, λ, lx, ly; quadrature=:rectangle)

Evaluate the first Rayleigh–Sommerfeld diffraction integral between parallel
planes using Shen and Wang's FFT direct integration (FFT-DI) method [1]. All
lengths in this matrix interface are in metres. The output is a
`Matrix{ComplexF64}` on the input grid, with pitches `lx/size(Uin, 1)` and
`ly/size(Uin, 2)`. The input is zero outside its finite window.

For `d > 0`, the full impulse response is
`h(x,y) = d*exp(2πim*r/λ)/(2π*r^2) * (1/r - 2πim/λ)`,
where `r = sqrt(x^2 + y^2 + d^2)`. Both the spherical phase and near-field
`1/r` term are retained. Zero-padded FFTs evaluate a linear convolution without
periodic wraparound. Sampling must resolve the input and the oscillatory kernel;
padding alone cannot correct undersampling, especially near `d = 0`.

`quadrature=:rectangle` assigns every input sample the area `Δx*Δy`, consistent
with the half-open physical windows used by [`MonoLightField2D`](@ref).
`quadrature=:simpson` uses the tensor product of Simpson's weights
`[1,4,2,…,4,1]/3` [1, Eqs. (17)–(20)]. It requires odd sample counts of at least
three on both axes and integrates only between the first and last sample
coordinates, spanning `(N-1)*Δ` rather than the full window width `N*Δ`.

`d` must be finite and nonnegative; `λ`, `lx` and `ly` must be finite and
positive. The input must be a nonempty one-based matrix. At zero distance,
return an independent complex copy, bypassing the singular kernel (but still
validating the grid and quadrature). Unlike scaled convolution [2], this method
preserves the input sampling and window size.

[1] F. Shen and A. Wang, "Fast-Fourier-transform based numerical integration
method for the Rayleigh-Sommerfeld diffraction formula," Applied Optics 45,
1102–1110 (2006), https://doi.org/10.1364/AO.45.001102, Eqs. (7)–(20).

[2] V. Nascov and P. C. Logofătu, "Fast computation algorithm for the
Rayleigh-Sommerfeld diffraction formula using a type of scaled convolution,"
Applied Optics 48, 4310–4319 (2009), https://doi.org/10.1364/AO.48.004310.
"""
function rayleighsommerfeld(Uin::AbstractMatrix{<:Number}, d::Real, λ::Real,
                           lx::Real, ly::Real; quadrature::Symbol=:rectangle)
    n, m = _propagation_grid(Uin, d, λ, lx, ly)
    d >= 0 || throw(ArgumentError("distance must be nonnegative"))
    quadrature in (:rectangle, :simpson) ||
        throw(ArgumentError("quadrature must be :rectangle or :simpson"))
    if quadrature === :simpson
        n >= 3 && m >= 3 && isodd(n) && isodd(m) ||
            throw(ArgumentError("Simpson quadrature requires odd sample counts of at least three on both axes"))
    end
    d == 0 && return Matrix{ComplexF64}(Uin)

    dx, dy = Float64(lx)/n, Float64(ly)/m
    # Smooth FFT lengths avoid expensive prime-size transforms. Each possible
    # source-to-target lag is stored in circular order; extra entries stay zero.
    np, mp = nextprod((2, 3, 5), 2n - 1), nextprod((2, 3, 5), 2m - 1)
    kernel = zeros(ComplexF64, np, mp)
    carrier = cispi(2*d/λ)
    for j in (1-m):(m-1), i in (1-n):(n-1)
        ρ = hypot(i*dx, j*dy)
        r = hypot(d, ρ)
        # r-d = ρ²/(r+d) avoids losing the transverse phase when d ≫ ρ.
        path_difference = ρ*(ρ/(r + d))
        phase = carrier * cispi(2*path_difference/λ)
        kernel[mod(i, np)+1, mod(j, mp)+1] =
            phase * (d/r)/r * (inv(2pi*r) - im/λ)
    end

    padded = zeros(ComplexF64, np, mp)
    input = Matrix{ComplexF64}(Uin)
    for j in 1:m, i in 1:n
        wx = quadrature === :simpson ? (i == 1 || i == n ? 1 : iseven(i) ? 4 : 2)/3 : 1.0
        wy = quadrature === :simpson ? (j == 1 || j == m ? 1 : iseven(j) ? 4 : 2)/3 : 1.0
        padded[i, j] = input[i, j]*wx*wy
    end
    fft!(padded)
    fft!(kernel)
    padded .*= kernel
    ifft!(padded)
    padded[1:n, 1:m] .* (dx*dy)
end

"""
    rayleighsommerfeld(Uin::MonoLightField2D, d::Unitful.Length;
                      quadrature=:rectangle) -> MonoLightField2D

Propagate a monochromatic light field with [`rayleighsommerfeld`](@ref),
preserving its wavelength, physical widths, sample counts and window centre.
Input and output windows may be translated together without changing the result.
"""
function rayleighsommerfeld(Uin::MonoLightField2D, d::Unitful.Length;
                           quadrature::Symbol=:rectangle)
    Uout = rayleighsommerfeld(Uin.data, auval(d), auval(Uin.wavelength),
                             auval(Uin.size[1]), auval(Uin.size[2]); quadrature=quadrature)
    _withdata(Uin, Uout)
end

"""
    rayleighsommerfeld(d; kwargs...)

Return a propagation functor with fixed distance and optional quadrature rule.
For example, `rayleighsommerfeld(1u"cm")(field)` propagates by 1 cm.
"""
rayleighsommerfeld(d; kwargs...) = x -> rayleighsommerfeld(x, d; kwargs...)
