"""
    angular_spectrum(Uin, d, λ, lx, ly)

calculate the propagation light field based on the angular spectrum.
"""
function angular_spectrum(Uin::AbstractArray{<:Number,2}, d::Real, λ::Real, lx::Real, ly::Real)
    tf = (u, v) -> 1-u^2-v^2>=0 ? exp(2im*pi*d/λ*sqrt(1-u^2-v^2)) : 0
    (n, m) = size(Uin)
    (u, v) = (-n/2+1:n/2,-m/2+1:m/2)./(lx,ly).*λ
    ifft(fft(Uin).*fftshift(tf.(u',v)))
end
