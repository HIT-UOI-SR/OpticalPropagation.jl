"""
    rayleighsommerfeld(Uin, d, λ, lx, ly; quadrature=:rectangle,
                      output_size=(lx, ly), output_samples=size(Uin), padding=1)

Evaluate the first Rayleigh–Sommerfeld diffraction integral between parallel
planes using FFT direct integration (FFT-DI) [1] or Nascov and Logofătu's scaled
convolution [2]. All lengths in this matrix interface are in metres. The output
is a `Matrix{ComplexF64}`. The input is zero outside its finite window.

`output_size` and `output_samples` independently specify the output window
widths and positive integer sample counts. Each defaults to the corresponding
input value. The pitches are `output_size ./ output_samples`. Both windows are
centred at zero, with half-open coordinates `x[i] = -Lx/2 + (i-1)*Lx/Nx`.
Equal input and output grids use ordinary FFT-DI; other grids use scaled
convolution with chirp-based scaled Fourier transforms [2, Eqs. (14)–(22) and
Appendix A], generalized to rectangular arrays and unequal sample counts.

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
validating the grids and quadrature). A different output grid at zero distance
is unsupported and raises an `ArgumentError`.

Scaled convolution uses a periodic Fourier interpolant of the sampled kernel,
not its exact values at fractional displacements. Its kernel pitches are the
larger input/output pitch on each axis. Both pitches must resolve the kernel;
this is especially demanding for small distances and coarse output sampling.
`padding` is a positive integer multiplying the kernel period in scaled mode.
Increasing it moves the periodic boundary away from the required displacements
but does not refine the kernel pitch. Check convergence with finer input and
output sampling. The unscaled FFT-DI path ignores `padding` after validation.

[1] F. Shen and A. Wang, "Fast-Fourier-transform based numerical integration
method for the Rayleigh-Sommerfeld diffraction formula," Applied Optics 45,
1102–1110 (2006), https://doi.org/10.1364/AO.45.001102, Eqs. (7)–(20).

[2] V. Nascov and P. C. Logofătu, "Fast computation algorithm for the
Rayleigh-Sommerfeld diffraction formula using a type of scaled convolution,"
Applied Optics 48, 4310–4319 (2009), https://doi.org/10.1364/AO.48.004310.
"""
function rayleighsommerfeld(Uin::AbstractMatrix{<:Number}, d::Real, λ::Real,
                           lx::Real, ly::Real; quadrature::Symbol=:rectangle,
                           output_size::Tuple{<:Real,<:Real}=(lx, ly),
                           output_samples::Tuple{<:Integer,<:Integer}=size(Uin),
                           padding::Integer=1)
    n, m = _propagation_grid(Uin, d, λ, lx, ly)
    d >= 0 || throw(ArgumentError("distance must be nonnegative"))
    quadrature in (:rectangle, :simpson) ||
        throw(ArgumentError("quadrature must be :rectangle or :simpson"))
    if quadrature === :simpson
        n >= 3 && m >= 3 && isodd(n) && isodd(m) ||
            throw(ArgumentError("Simpson quadrature requires odd sample counts of at least three on both axes"))
    end
    all(x -> isfinite(x) && x > 0, output_size) ||
        throw(ArgumentError("output widths must be finite and positive"))
    all(>(0), output_samples) || throw(ArgumentError("output sample counts must be positive"))
    padding > 0 || throw(ArgumentError("padding must be positive"))
    samegrid = output_size == (lx, ly) && output_samples == (n, m)
    if d == 0
        samegrid || throw(ArgumentError("a different output grid is unsupported at zero distance"))
        return Matrix{ComplexF64}(Uin)
    end

    input = Matrix{ComplexF64}(Uin)
    if quadrature === :simpson
        for j in 1:m, i in 1:n
            wx = (i == 1 || i == n ? 1 : iseven(i) ? 4 : 2)/3
            wy = (j == 1 || j == m ? 1 : iseven(j) ? 4 : 2)/3
            input[i, j] *= wx*wy
        end
    end
    if !samegrid
        return _rs_scaled_convolution(input, d, λ, (lx, ly), output_size, output_samples, padding)
    end

    dx, dy = Float64(lx)/n, Float64(ly)/m
    # Smooth FFT lengths avoid expensive prime-size transforms. Each possible
    # source-to-target lag is stored in circular order; extra entries stay zero.
    np, mp = nextprod((2, 3, 5), 2n - 1), nextprod((2, 3, 5), 2m - 1)
    kernel = zeros(ComplexF64, np, mp)
    carrier = cispi(2*d/λ)
    for j in (1-m):(m-1), i in (1-n):(n-1)
        kernel[mod(i, np)+1, mod(j, mp)+1] =
            _rs_kernel(i*dx, j*dy, d, λ, carrier)
    end

    padded = zeros(ComplexF64, np, mp)
    padded[1:n, 1:m] .= input
    fft!(padded)
    fft!(kernel)
    padded .*= kernel
    ifft!(padded)
    padded[1:n, 1:m] .* (dx*dy)
end

function _rs_kernel(x, y, d, λ, carrier)
    ρ = hypot(x, y)
    r = hypot(d, ρ)
    # r-d = ρ²/(r+d) avoids losing the transverse phase when d ≫ ρ.
    path_difference = ρ*(ρ/(r + d))
    phase = carrier * cispi(2*path_difference/λ)
    phase * (d/r)/r * (inv(2pi*r) - im/λ)
end

# Along one axis, evaluate sum_i A[i] exp(2πim*β*(i+a)*(j+b)),
# where i,j are zero-based and a,b can be fractional coordinate origins.
# Bluestein's identity 2ij = i²+j²-(j-i)² turns the rectangular transform
# into a linear chirp convolution, as in Nascov and Logofătu's Appendix A.
function _rs_scaled_dft(A, dim, nout, β, a, b)
    nin = size(A, dim)
    nf = nextprod((2, 3, 5), nin + nout - 1)
    shape = dim == 1 ? (nf, size(A, 2)) : (size(A, 1), nf)
    reshape_axis(v) = reshape(v, dim == 1 ? (length(v), 1) : (1, length(v)))
    padded = zeros(ComplexF64, shape)
    src = dim == 1 ? view(padded, 1:nin, :) : view(padded, :, 1:nin)
    pre = [cispi(β*i^2 + 2β*b*i) for i in 0:nin-1]
    src .= A .* reshape_axis(pre)
    chirp = zeros(ComplexF64, nf)
    for k in (1-nin):(nout-1)
        chirp[mod(k, nf)+1] = cispi(-β*k^2)
    end
    fft!(chirp)
    fft!(padded, (dim,))
    padded .*= reshape_axis(chirp)
    ifft!(padded, (dim,))
    post = [cispi(β*j^2 + 2β*a*(j+b)) for j in 0:nout-1]
    dst = dim == 1 ? view(padded, 1:nout, :) : view(padded, :, 1:nout)
    dst .* reshape_axis(post)
end

function _rs_scaled_convolution(input, d, λ, input_size, output_size, output_samples, padding)
    ni = size(input)
    no = output_samples
    di = Float64.(input_size) ./ ni
    do_ = Float64.(output_size) ./ no
    δ = max.(di, do_)
    # A centred even period encloses every input-to-output displacement.
    # This is Eq. (17) extended to unequal counts and half-open odd grids.
    periods = ntuple(2) do axis
        halfwidth = (input_size[axis]/2 + output_size[axis]/2)/δ[axis]
        2*padding*nextprod((2, 3, 5), ceil(Int, halfwidth) + 1)
    end
    nx, ny = periods
    kernel = Matrix{ComplexF64}(undef, nx, ny)
    carrier = cispi(2*d/λ)
    for j in 1:ny, i in 1:nx
        x = (i-1 < nx÷2 ? i-1 : i-1-nx)*δ[1]
        y = (j-1 < ny÷2 ? j-1 : j-1-ny)*δ[2]
        kernel[i, j] = _rs_kernel(x, y, d, λ, carrier)
    end
    fft!(kernel)
    spectrum = fftshift(kernel)
    # Physical coordinates use i-N/2, including half-integer origins for odd N.
    transformed = _rs_scaled_dft(input, 1, nx, -di[1]/(δ[1]*nx), -ni[1]/2, -nx/2)
    transformed = _rs_scaled_dft(transformed, 2, ny, -di[2]/(δ[2]*ny), -ni[2]/2, -ny/2)
    transformed .*= spectrum
    output = _rs_scaled_dft(transformed, 2, no[2], do_[2]/(δ[2]*ny), -ny/2, -no[2]/2)
    output = _rs_scaled_dft(output, 1, no[1], do_[1]/(δ[1]*nx), -nx/2, -no[1]/2)
    output .*= (di[1]*di[2])/(nx*ny)
end

"""
    rayleighsommerfeld(Uin::MonoLightField2D, d::Unitful.Length;
                      quadrature=:rectangle, output_size=physicalsize(Uin),
                      output_samples=size(Uin), padding=1) -> MonoLightField2D

Propagate a monochromatic light field with [`rayleighsommerfeld`](@ref),
preserving its wavelength and window centre. `output_size` is a pair of lengths
with units; `output_samples` sets the output sample counts independently. By
default the physical widths and counts are preserved. Input and output windows
may be translated together without changing the result.
"""
function rayleighsommerfeld(Uin::MonoLightField2D, d::Unitful.Length;
                           quadrature::Symbol=:rectangle,
                           output_size::Tuple{<:Unitful.Length,<:Unitful.Length}=physicalsize(Uin),
                           output_samples::Tuple{<:Integer,<:Integer}=size(Uin),
                           padding::Integer=1)
    Uout = rayleighsommerfeld(Uin.data, auval(d), auval(Uin.wavelength),
                             auval(Uin.size[1]), auval(Uin.size[2]); quadrature=quadrature,
                             output_size=auval.(output_size), output_samples=output_samples,
                             padding=padding)
    _withdata(Uin, Uout; size=output_size)
end

"""
    rayleighsommerfeld(d; kwargs...)

Return a propagation functor with fixed distance, quadrature and output grid.
For example, `rayleighsommerfeld(1u"cm")(field)` propagates by 1 cm.
"""
rayleighsommerfeld(d; kwargs...) = x -> rayleighsommerfeld(x, d; kwargs...)
