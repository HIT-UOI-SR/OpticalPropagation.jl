"""
    fresnel2(Uin, d, λ, lx, ly)

calculate the propagation light field based on the Fresnel diffraction with double Fourier transform.

The input and output share the same uniformly sampled periodic window. The
transfer function uses FFT-order frequencies for any odd/even sample counts.
`d` may be positive or negative and must be finite. Zero distance returns an
independent complex copy. Wavelength and window widths must be finite and positive.

## Arguments

- `Uin::AbstractArray{<:Number,2}`: Complex array of input complex amplitude.
- `d::Real`: Distance to propagate in metres.
- `λ::Real`: Wavelength of light to propagate.
- `lx::Real`: The physical size of the input data along the x-axis.
- `ly::Real`: The physical size of the input data along the y-axis.

## Returns

- `::Array{<:Number,2}`: Complex amplitude data after propagation.
"""
function fresnel2(Uin::AbstractArray{<:Number,2}, d::Real, λ::Real, lx::Real, ly::Real)
    n, m = _propagation_grid(Uin, d, λ, lx, ly)
    spectrum = Matrix{ComplexF64}(Uin)
    d == 0 && return spectrum
    fft!(spectrum)
    fx, fy = fftfreq(n, n/lx), fftfreq(m, m/ly)
    hx = [cispi(-λ*d*f^2) for f in fx]
    hy = [cispi(-λ*d*f^2) for f in fy]
    carrier = cispi(2*d/λ)
    for j in 1:m, i in 1:n
        spectrum[i, j] *= carrier * hx[i] * hy[j]
    end
    ifft!(spectrum)
end

"""
    fresnel2(Uin::MonoLightField2D, d::Unitful.Length) -> MonoLightField2D

calculate the propagation light field for [`MonoLightField2D`](@ref) based on the Fresnel diffraction with double Fourier transform.
"""
function fresnel2(Uin::MonoLightField2D, d::Unitful.Length)
    Uout = fresnel2(Uin.data, auval(d), auval(Uin.wavelength), auval(Uin.size[1]), auval(Uin.size[2]))
    _withdata(Uin, Uout)
end

"""
    fresnel2(d)

represent a [`fresnel2`](@ref) propagation functor with the specified distance.
"""
fresnel2(d) = x->fresnel2(x, d)

"""
    fresnel1(Uin, d, λ, lx, ly)

calculate the propagation light field based on the Fresnel diffraction with single Fourier transform.

Input coordinates are `x[i] = -lx/2 + (i-1)*lx/n`, and similarly for y.
The output follows the same half-open convention with positive widths
`lxo = λ*abs(d)/(lx/n)` and `lyo = λ*abs(d)/(ly/m)`. Both quadratic phases and
the phase from the coordinate origins are included, for odd/even and singleton
dimensions. A negative distance uses an unnormalised backward FFT and retains
increasing output coordinates. Zero distance returns a copy with the input widths.
`d` must be finite; wavelength and input widths must be finite and positive.

## Arguments

- `Uin::AbstractArray{<:Number,2}`: Complex array of input complex amplitude.
- `d::Real`: Distance to propagate in metres.
- `λ::Real`: Wavelength of light to propagate.
- `lx::Real`: The physical size of the input data along the x-axis.
- `ly::Real`: The physical size of the input data along the y-axis.

## Returns

    (Uout, (lxo, lyo))

- `Uout::Array{<:Number,2}`: Complex amplitude data after propagation.
- `lxo::Real`: The physical size of the diffraction data along the x-axis.
- `lyo::Real`: The physical size of the diffraction data along the y-axis.
"""
function fresnel1(Uin::AbstractArray{<:Number,2}, d::Real, λ::Real, lx::Real, ly::Real)
    n, m = _propagation_grid(Uin, d, λ, lx, ly)
    d == 0 && return (Matrix{ComplexF64}(Uin), (float(lx), float(ly)))
    dx = lx/n
    dy = ly/m
    lxo = λ*abs(d)/dx
    lyo = λ*abs(d)/dy
    chirp_in_x, chirp_out_x = _fresnel1_chirps(n, dx, lxo/n, λ*d)
    chirp_in_y, chirp_out_y = _fresnel1_chirps(m, dy, lyo/m, λ*d)
    amplitude = Matrix{ComplexF64}(undef, n, m)
    for j in 1:m, i in 1:n
        amplitude[i, j] = Uin[i, j] * chirp_in_x[i] * chirp_in_y[j]
    end
    d > 0 ? fft!(amplitude) : bfft!(amplitude)
    # For zero-based indices a,b, the centred cross term is
    # exp(-2πi*sign(d)*(a-n/2)*(b-n/2)/n). Its expansion gives the
    # checkerboard factors in the chirps and this constant origin phase.
    origin_phase = cispi(-sign(d)*(n + m)/2)
    scale = cispi(2*d/λ) * origin_phase * dx*dy/(im*λ*d)
    for j in 1:m, i in 1:n
        amplitude[i, j] *= scale * chirp_out_x[i] * chirp_out_y[j]
    end
    (amplitude, (lxo, lyo))
end

function _fresnel1_chirps(n, pitch_in, pitch_out, λd)
    incoming = Vector{ComplexF64}(undef, n)
    outgoing = similar(incoming)
    for i in 1:n
        offset = i - 1 - n/2
        modulation = isodd(i - 1) ? -1 : 1
        incoming[i] = modulation * cispi((offset*pitch_in)^2/λd)
        outgoing[i] = modulation * cispi((offset*pitch_out)^2/λd)
    end
    incoming, outgoing
end

"""
    fresnel1(Uin::MonoLightField2D, d::Unitful.Length) -> MonoLightField2D

calculate the propagation light field for [`MonoLightField2D`](@ref) based on the Fresnel diffraction with single Fourier transform.
"""
function fresnel1(Uin::MonoLightField2D, d::Unitful.Length)
    (Uout, Lout) = fresnel1(Uin.data, auval(d), auval(Uin.wavelength), auval(Uin.size[1]), auval(Uin.size[2]))
    _withdata(Uin, Uout; size=Lout .* u"m")
end

"""
    fresnel1(d)

represent a [`fresnel1`](@ref) propagation functor with the specified distance.
"""
fresnel1(d) = x->fresnel1(x, d)
