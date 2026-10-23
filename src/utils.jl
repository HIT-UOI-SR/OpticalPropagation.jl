function auval(x::Unitful.Length) # get absolute numerical value of quantity
    unit = u"m"
    convert(Float64, x / unit)
end

# Shared domain for the planar FFT propagators. Distances may be signed, while
# wavelengths and physical window widths must be positive.
function _propagation_grid(Uin, d::Real, λ::Real, lx::Real, ly::Real)
    Base.require_one_based_indexing(Uin)
    isfinite(d) || throw(ArgumentError("distance must be finite"))
    isfinite(λ) && λ > 0 || throw(ArgumentError("wavelength must be finite and positive"))
    isfinite(lx) && lx > 0 || throw(ArgumentError("lx must be finite and positive"))
    isfinite(ly) && ly > 0 || throw(ArgumentError("ly must be finite and positive"))
    n, m = size(Uin)
    n > 0 && m > 0 || throw(ArgumentError("the input field must be nonempty"))
    n, m
end
