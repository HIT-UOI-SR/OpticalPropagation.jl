"""
    fresnel_diffraction_double(data, d, λ, lx, ly)

calculate the propagation light field based on the Fresnel diffraction with double Fourier transform.
"""
function fresnel_diffraction_double(data::AbstractArray{<:Number,2}, d::Real, λ::Real, lx::Real, ly::Real)
    tf = (u, v) -> exp(2im*pi*d/λ)*exp(-im*pi*d/λ*(u^2+v^2))
    (n, m) = size(data)
    (u, v) = (-n/2+1:n/2,-m/2+1:m/2)./(lx,ly).*λ
    ifft(fft(data).*fftshift(tf.(u',v)))
end

"""
    fresnel_diffraction_single(data, d, λ, lx, ly)

calculate the propagation light field based on the Fresnel diffraction with single Fourier transform.
"""
function fresnel_diffraction_single(data::AbstractArray{<:Number,2}, d::Real, λ::Real, lx::Real, ly::Real)
    qp = (x, y) -> exp(im*pi/d/λ*(x^2+y^2))
    (n, m) = size(data)
    (dx, dy) = (lx, ly)./(n, m)
    (lxo, lyo) = λ*d./(dx, dy)
    (u, v) = (-n/2+1:n/2,-m/2+1:m/2)./(lx,ly).*λ
    amp = fft(data.*qp.(u', v))
    (uo, vo) = (-n/2+1:n/2,-m/2+1:m/2)./(lxo,lyo).*λ
    patt = amp.*qp.(u', v).*exp(2im*pi*d/λ).*dx.*dy./(im*λ*d)
    (patt, (lxo, lyo))
end
