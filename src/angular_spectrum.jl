"""
    angular_spectrum(data, d, λ, lx, ly)

calculate the propagation light field based on the angular spectrum.
"""
function angular_spectrum(data::AbstractArray{<:Number,2}, d::Real, λ::Real, lx::Real, ly::Real)
    tf = (u, v) -> 1-u^2-v^2>=0 ? exp(2im*pi*d/λ*sqrt(1-u^2-v^2)) : 0
    (n, m) = size(data)
    (u, v) = (-n/2+1:n/2,-m/2+1:m/2)./(lx,ly).*λ
    ifft(fft(data).*fftshift(tf.(u',v)))
end
