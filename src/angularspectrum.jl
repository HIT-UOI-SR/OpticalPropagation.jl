"""
    angularspectrum(Uin, d, λ, lx, ly)

calculate the propagation light field based on the angular spectrum.

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
    (n, m) = size(Uin)
    ua = (-n/2+1:n/2)/lx*λ
    va = (-m/2+1:m/2)/ly*λ
    ifft(fft(Uin).*fftshift([1.0-u^2-v^2>=0.0 ? exp(2im*pi*d/λ*sqrt(1.0-u^2-v^2)) : 0.0+0.0im for v in va, u in ua]))
end


"""
    angularspectrum(Uin::MonoLightField2D, d::Unitful.Length) -> MonoLightField2D

calculate the propagation light field for [`MonoLightField2D`](@ref) based on the angular spectrum.
"""
function angularspectrum(Uin::MonoLightField2D, d::Unitful.Length)
    Uout = angularspectrum(Uin.data, uval(d), uval(Uin.wavelength), uval(Uin.size[1]), uval(Uin.size[2]))
    MonoLightField2D(Uin, data=Uout)
end
