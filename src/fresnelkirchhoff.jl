"""
    fresnelkirchhoff(Uin, d, λ, lx, ly; radius=Inf)

Propagate a sampled scalar field with a Fresnel–Kirchhoff zone-mask convolution,
following Hudson's spherical-reference-surface approximation [1]. All lengths in
this array interface are in metres. The result is a `Matrix{ComplexF64}` with the
same number of samples as `Uin`; input samples outside the array are zero.

`d` is a finite, nonnegative axial distance; `λ`, `lx`, and `ly` are finite and
positive. `lx` and `ly` are the input widths along the first and second array
dimensions, with sample pitches `lx/size(Uin, 1)` and `ly/size(Uin, 2)`.
Zero distance returns an independent complex copy of the input.

By default, the input and output are parallel planes with the same physical size,
and normally incident illumination is assumed. The impulse response is
`exp(2πim*s/λ) * (1 + d/s) / (2im*λ*s)`, where `s = sqrt(d^2 + x^2 + y^2)`.
The FFTs are padded to compute a linear convolution without periodic wraparound.

A finite, nonzero `radius` specifies a signed input reference-sphere radius:
positive for a diverging wave and negative for a converging wave. The output
sphere is concentric, with signed radius `radius + d`. Its projected physical
size is `(1 + d/radius) .* (lx, ly)`. Propagation to or through the common centre
is unsupported. The data give the actual complex amplitude on the input cap,
including its reference-wave amplitude and phase. The returned data similarly
give the actual amplitude on the output cap.

The grid is projected perpendicularly from the tangent plane onto the sphere,
with coordinates `x[i] = -lx/2 + (i-1)*lx/n` and similarly for `y`, as in the
light-field plotting recipe. The surface-area Jacobian is included. Both the
input and the entire displacement grid must fit inside a hemisphere.
Translation on this grid approximates rotation on the sphere, so use caps small
compared with `abs(radius)`. The oscillatory zone mask must also be adequately
sampled; zero padding does not correct undersampling or aperture truncation.

[1] J. A. Hudson, "Fresnel-Kirchhoff diffraction in optical systems: an approximate
computational algorithm," Applied Optics 23, 2292–2295 (1984),
https://doi.org/10.1364/AO.23.002292, Eqs. (1)–(5).
"""
function fresnelkirchhoff(Uin::AbstractMatrix{<:Number}, d::Real, λ::Real,
                         lx::Real, ly::Real; radius::Real=Inf)
    isfinite(d) && d >= 0 || throw(ArgumentError("distance must be finite and nonnegative"))
    isfinite(λ) && λ > 0 || throw(ArgumentError("wavelength must be finite and positive"))
    isfinite(lx) && lx > 0 || throw(ArgumentError("lx must be finite and positive"))
    isfinite(ly) && ly > 0 || throw(ArgumentError("ly must be finite and positive"))
    radius == Inf || (isfinite(radius) && radius != 0) ||
        throw(ArgumentError("radius must be finite and nonzero, or +Inf for a plane"))
    1 + d/radius > 0 || throw(ArgumentError("propagation to or through the sphere centre is unsupported"))
    n, m = size(Uin)
    n > 0 && m > 0 || throw(ArgumentError("the input field must be nonempty"))
    dx, dy = Float64(lx)/n, Float64(ly)/m
    hypot(max(n/2, n - 1)*dx, max(m/2, m - 1)*dy) < abs(radius) ||
        throw(ArgumentError("the input and displacement grids must fit strictly inside a reference hemisphere"))
    d == 0 && return Matrix{ComplexF64}(Uin)

    # FFT order: zero displacement first, then positive and negative offsets.
    # 2n-1 by 2m-1 samples hold every source-to-target displacement exactly once.
    np, mp = 2n - 1, 2m - 1
    kernel = Matrix{ComplexF64}(undef, np, mp)
    for j in 1:mp, i in 1:np
        x = (i <= n ? i - 1 : i - 1 - np)*dx
        y = (j <= m ? j - 1 : j - 1 - mp)*dy
        ρ = hypot(x, y)
        c = sqrt(1 - (ρ/radius)^2)
        # Stable sphere sag, including the planar limit and signed curvature.
        sag = (ρ/radius)*ρ/(1 + c)
        s = hypot(ρ, d + sag)
        cosχ = (d*c - sag)/s
        kernel[i, j] = cispi(2s/λ)*(1 + cosχ)/(2im*λ*s)
    end

    padded = zeros(ComplexF64, np, mp)
    input = Matrix{ComplexF64}(Uin)
    for j in 1:m, i in 1:n
        x = (i - 1 - n/2)*dx
        y = (j - 1 - m/2)*dy
        # dA = dx*dy/cos(θ) for orthogonal projection onto the input cap.
        c = sqrt(1 - (hypot(x, y)/radius)^2)
        padded[i, j] = input[i, j]/c
    end
    fft!(padded)
    fft!(kernel)
    padded .*= kernel
    ifft!(padded)
    padded[1:n, 1:m] .* (dx*dy)
end

"""
    fresnelkirchhoff(Uin::MonoLightField2D, d::Unitful.Length;
                     radius::Unitful.Length=Inf*u"m") -> MonoLightField2D

Propagate a monochromatic light field with [`fresnelkirchhoff`](@ref). A finite
signed `radius` selects a spherical reference surface and scales the output
physical size by `1 + d/radius`. The wavelength is preserved.
The planar mode preserves the window centre. Finite-radius propagation currently
requires a zero-centred window; support for off-axis spherical caps is deferred.
"""
function fresnelkirchhoff(Uin::MonoLightField2D, d::Unitful.Length;
                         radius::Unitful.Length=Inf*u"m")
    if isfinite(radius) && any(!iszero, gridcenter(Uin))
        throw(ArgumentError("finite-radius propagation currently requires a zero-centred window"))
    end
    Uout = fresnelkirchhoff(Uin.data, auval(d), auval(Uin.wavelength),
                          auval(Uin.size[1]), auval(Uin.size[2]); radius=auval(radius))
    scale = 1 + auval(d)/auval(radius)
    _withdata(Uin, Uout; size=scale .* Uin.size, center=scale .* gridcenter(Uin))
end

"""
    fresnelkirchhoff(d; kwargs...)

Return a propagation functor with fixed distance and optional reference radius.
For example, `fresnelkirchhoff(1u"cm")(field)` propagates a light field by 1 cm.
"""
fresnelkirchhoff(d; kwargs...) = x -> fresnelkirchhoff(x, d; kwargs...)
