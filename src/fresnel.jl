"""
    fresnel2(Uin, d, λ, lx, ly)

calculate the propagation light field based on the Fresnel diffraction with double Fourier transform.

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
    tf = (u, v) -> exp(2im*pi*d/λ)*exp(-im*pi*d/λ*(u^2+v^2))
    (n, m) = size(Uin)
    u = (-n/2+1:n/2)/lx*λ
    v = (-m/2+1:m/2)/ly*λ
    ifft(fft(Uin).*fftshift(tf.(u',v)))
end

"""
    fresnel2(Uin::MonoLightField2D, d::Unitful.Length) -> MonoLightField2D

calculate the propagation light field for [`MonoLightField2D`](@ref) based on the Fresnel diffraction with double Fourier transform.
"""
function fresnel2(Uin::MonoLightField2D, d::Unitful.Length)
    Uout = fresnel2(Uin.data, uval(d), uval(Uin.wavelength), uval(Uin.size[1]), uval(Uin.size[2]))
    MonoLightField2D(Uin, data=Uout)
end

"""
    fresnel1(Uin, d, λ, lx, ly)

calculate the propagation light field based on the Fresnel diffraction with single Fourier transform.

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
    qp = (x, y) -> exp(im*pi/d/λ*(x^2+y^2))
    (n, m) = size(Uin)
    dx = lx/n
    dy = ly/m
    lxo = λ*d/dx
    lyo = λ*d/dy
    u = (-n/2+1:n/2)/lx*λ
    v = (-m/2+1:m/2)/ly*λ
    amp = fft(Uin.*qp.(u', v))
    u *= lx/lxo
    v *= ly/lyo
    Uout = amp.*qp.(u', v).*exp(2im*pi*d/λ).*dx.*dy./(im*λ*d)
    (Uout, (lxo, lyo))
end

"""
    fresnel1(Uin::MonoLightField2D, d::Real) -> MonoLightField2D

calculate the propagation light field for [`MonoLightField2D`](@ref) based on the Fresnel diffraction with single Fourier transform.
"""
function fresnel1(Uin::MonoLightField2D, d::Unitful.Length)
    (Uout, Lout) = fresnel1(Uin.data, uval(d), uval(Uin.wavelength), uval(Uin.size[1]), uval(Uin.size[2]))
    MonoLightField2D(Uin, data=Uout, size=Lout)
end
