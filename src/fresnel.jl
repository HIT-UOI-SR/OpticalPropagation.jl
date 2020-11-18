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
    (n, m) = size(Uin)
    ua = (-n/2:n/2-1)/lx*λ
    va = (-m/2:m/2-1)/ly*λ
    ifft(fft(Uin).*ifftshift([exp(2im*pi*d/λ)*exp(-im*pi*d/λ*(u^2+v^2)) for u in ua, v in va]))
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
    (n, m) = size(Uin)
    dx = lx/n
    dy = ly/m
    lxo = λ*d/dx
    lyo = λ*d/dy
    xa = (-n/2:n/2-1)/lx*λ
    ya = (-m/2:m/2-1)/ly*λ
    amp = fftshift(fft(Uin.*[exp(im*pi/d/λ*(x^2+y^2)) for x in xa, y in ya]))
    xa *= lx/lxo
    ya *= ly/lyo
    Uout = amp.*[exp(im*pi/d/λ*(x^2+y^2)) for x in xa, y in ya].*exp(2im*pi*d/λ).*dx.*dy./(im*λ*d)
    (Uout, (lxo, lyo))
end

"""
    fresnel1(Uin::MonoLightField2D, d::Unitful.Length) -> MonoLightField2D

calculate the propagation light field for [`MonoLightField2D`](@ref) based on the Fresnel diffraction with single Fourier transform.
"""
function fresnel1(Uin::MonoLightField2D, d::Unitful.Length)
    (Uout, Lout) = fresnel1(Uin.data, uval(d), uval(Uin.wavelength), uval(Uin.size[1]), uval(Uin.size[2]))
    MonoLightField2D(Uin, data=Uout, size=Lout)
end
