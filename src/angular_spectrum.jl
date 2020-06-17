"""
    angular_spectrum(Uin, d, λ, lx, ly)

calculate the propagation light field based on the angular spectrum.

## Parameters

- `Uin::AbstractArray{<:Number,2}`: Complex array of input complex amplitude.
- `d::Real`: Distance to propagate in metres.
- `λ::Real`: Wavelength of light to propagate.
- `lx::Real`: The physical size of the input data along the x-axis.
- `ly::Real`: The physical size of the input data along the y-axis.

## Returns

- `::Array{<:Number,2}`: Complex amplitude data after propagation.
"""
function angular_spectrum(Uin::AbstractArray{<:Number,2}, d::Real, λ::Real, lx::Real, ly::Real)
    tf = (u, v) -> 1-u^2-v^2>=0 ? exp(2im*pi*d/λ*sqrt(1-u^2-v^2)) : 0
    (n, m) = size(Uin)
    (u, v) = (-n/2+1:n/2,-m/2+1:m/2)./(lx,ly).*λ
    ifft(fft(Uin).*fftshift(tf.(u',v)))
end
