"""
    AbstractLightField

Abstract supertype of all 2-dimensional light field types.
"""
abstract type AbstractLightField2D end

+(x::AbstractLightField2D,y::AbstractLightField2D) = +(promote(x,y)...)
-(x::AbstractLightField2D,y::AbstractLightField2D) = -(promote(x,y)...)
Base.isapprox(x::AbstractLightField2D,y::AbstractLightField2D; kwargs...) = isapprox(promote(x,y)...,kwargs...)

const LengthType=typeof(1.0*u"m");


"""
    MonoLightField2D  <: AbstractLightField2D

`MonoLightField2D` is used to represent a 2-dimensional monochromatic light field.

# Examples

```jldoctest
julia> using Unitful

julia> MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",physical_size=(1u"mm",1u"mm"))
MonoLightField2D:
    wavelength: 6.328e-7 m
    physical size: (0.001 m, 0.001 m)
    distribution data:
2×2 Array{Complex{Float64},2}:
 1.0+0.0im  2.0+0.0im
 3.0+0.0im  4.0+0.0im
```
"""
struct MonoLightField2D  <: AbstractLightField2D
    distribution_data::Array{ComplexF64,2}
    wavelength::LengthType
    physical_size::Tuple{LengthType,LengthType}

    MonoLightField2D(data::AbstractMatrix; wavelength::Unitful.Length, physical_size::Tuple{<: Unitful.Length,<: Unitful.Length}) =
        new(ComplexF64.(data), float(wavelength), float.(physical_size))
end

function MonoLightField2D(x::MonoLightField2D; kwargs...)
    dist_data=get(kwargs,:distribution_data,x.distribution_data)
    wavelen=get(kwargs,:wavelength,x.wavelength)
    phy_size=get(kwargs,:physical_size,x.physical_size)
    MonoLightField2D(dist_data,wavelength=wavelen,physical_size=phy_size)
end

function Base.show(io::IO, x::MonoLightField2D)
    print(io, "MonoLightField2D(")
    show(io, x.distribution_data)
    print(io, ", wavelength=")
    show(io, x.wavelength)
    print(io, ", physical_size=")
    show(io, x.physical_size)
    print(io, ")")
end

function Base.show(io::IO, mime::MIME"text/plain", x::MonoLightField2D)
    print(io,"MonoLightField2D:")
    print(io,"\n    wavelength: ")
    show(io, mime, x.wavelength)
    print(io, "\n    physical size:")
    show(io, mime, x.physical_size)
    print(io,"\n    distribution data:\n")
    show(io, mime, x.distribution_data)
end

function +(a::MonoLightField2D,b::MonoLightField2D)
    if ! (a.wavelength≈b.wavelength)
        @error "Wavelengths $(a.wavelength) and $(b.wavelength) are not equal."
    elseif ! all(a.physical_size .≈ b.physical_size)
        @error "Physical size $(a.physical_size) and $(b.physical_size) are not equal."
    else
        MonoLightField2D(a.distribution_data .+ b.distribution_data,
            wavelength=a.wavelength,
            physical_size=a.physical_size)
    end
end

function -(a::MonoLightField2D,b::MonoLightField2D)
    if ! (a.wavelength≈b.wavelength)
        @error "Wavelengths $(a.wavelength) and $(b.wavelength) are not equal."
    elseif ! all(a.physical_size .≈ b.physical_size)
        @error "Physical size $(a.physical_size) and $(b.physical_size) are not equal."
    else
        MonoLightField2D(a.distribution_data .- b.distribution_data,
            wavelength=a.wavelength,
            physical_size=a.physical_size)
    end
end

*(r::Number, x::MonoLightField2D) =
    MonoLightField2D(r .* x.distribution_data,
        wavelength = x.wavelength,
        physical_size = x.physical_size)

*(x::MonoLightField2D, r::Number) = r*x

/(x::MonoLightField2D, r::Number) = (1/r)*x

\(r::Number, x::MonoLightField2D) = (1/r)*x

function Base.isapprox(a::MonoLightField2D,b::MonoLightField2D; kwargs...)
    isapprox(a.wavelength,b.wavelength,kwargs...) &&
      all(isapprox.(a.physical_size,b.physical_size,kwargs...)) &&
      isapprox(a.distribution_data,b.distribution_data,kwargs...)
end
