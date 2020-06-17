"""
    fresnel_diffraction_double(Uin, d, λ, lx, ly)

calculate the propagation light field based on the Fresnel diffraction with double Fourier transform.
"""
function fresnel_diffraction_double(Uin::AbstractArray{<:Number,2}, d::Real, λ::Real, lx::Real, ly::Real)
    tf = (u, v) -> exp(2im*pi*d/λ)*exp(-im*pi*d/λ*(u^2+v^2))
    (n, m) = size(Uin)
    (u, v) = (-n/2+1:n/2,-m/2+1:m/2)./(lx,ly).*λ
    ifft(fft(Uin).*fftshift(tf.(u',v)))
end

"""
    fresnel_diffraction_single(Uin, d, λ, lx, ly)

calculate the propagation light field based on the Fresnel diffraction with single Fourier transform.
"""
function fresnel_diffraction_single(Uin::AbstractArray{<:Number,2}, d::Real, λ::Real, lx::Real, ly::Real)
    qp = (x, y) -> exp(im*pi/d/λ*(x^2+y^2))
    (n, m) = size(Uin)
    (dx, dy) = (lx, ly)./(n, m)
    (lxo, lyo) = λ*d./(dx, dy)
    (u, v) = (-n/2+1:n/2,-m/2+1:m/2)./(lx,ly).*λ
    amp = fft(Uin.*qp.(u', v))
    (uo, vo) = (-n/2+1:n/2,-m/2+1:m/2)./(lxo,lyo).*λ
    Uout = amp.*qp.(u', v).*exp(2im*pi*d/λ).*dx.*dy./(im*λ*d)
    (Uout, (lxo, lyo))
end
