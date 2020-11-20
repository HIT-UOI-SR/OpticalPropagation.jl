"""
    AbstractLightField

Abstract supertype of all 2-dimensional light field types.
"""
abstract type AbstractLightField2D end

+(x::AbstractLightField2D,y::AbstractLightField2D) = +(promote(x,y)...)
-(x::AbstractLightField2D,y::AbstractLightField2D) = -(promote(x,y)...)
Base.isapprox(x::AbstractLightField2D,y::AbstractLightField2D; kwargs...) = isapprox(promote(x,y)...,kwargs...)

const LengthType=typeof(1.0u"m");


"""
    MonoLightField2D  <: AbstractLightField2D

`MonoLightField2D` is used to represent a 2-dimensional monochromatic light field.

# Examples

```jldoctest; setup = :(using OpticalPropagation)
julia> using Unitful

julia> MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
MonoLightField2D:
    wavelength: 6.328e-7 m
    physical size: (0.001 m, 0.001 m)
    complex amplitude data:
2×2 Array{Complex{Float64},2}:
 1.0+0.0im  2.0+0.0im
 3.0+0.0im  4.0+0.0im
```
"""
struct MonoLightField2D  <: AbstractLightField2D
    data::Array{ComplexF64,2}
    wavelength::LengthType
    size::Tuple{LengthType,LengthType}

    MonoLightField2D(data::AbstractMatrix; wavelength::Unitful.Length, size::Tuple{<: Unitful.Length,<: Unitful.Length}) =
        new(ComplexF64.(data), float(wavelength), float.(size))
end

function MonoLightField2D(x::MonoLightField2D; kwargs...)
    dist_data=get(kwargs,:data,x.data)
    wavelen=get(kwargs,:wavelength,x.wavelength)
    phy_size=get(kwargs,:size,x.size)
    MonoLightField2D(dist_data,wavelength=wavelen,size=phy_size)
end

function Base.show(io::IO, x::MonoLightField2D)
    print(io, "MonoLightField2D(")
    show(io, x.data)
    print(io, ", wavelength=")
    show(io, x.wavelength)
    print(io, ", size=")
    show(io, x.size)
    print(io, ")")
end

function Base.show(io::IO, mime::MIME"text/plain", x::MonoLightField2D)
    print(io,"MonoLightField2D:")
    print(io,"\n    wavelength: ")
    show(io, mime, x.wavelength)
    print(io, "\n    physical size: ")
    show(io, mime, x.size)
    print(io,"\n    complex amplitude data:\n")
    show(io, mime, x.data)
end

function +(a::MonoLightField2D,b::MonoLightField2D)
    if ! (a.wavelength≈b.wavelength)
        throw(DimensionMismatch("Wavelengths $(a.wavelength) and $(b.wavelength) are not equal."))
    elseif ! all(a.size .≈ b.size)
        throw(DimensionMismatch("Physical size $(a.size) and $(b.size) are not equal."))
    else
        MonoLightField2D(a.data .+ b.data,
            wavelength=a.wavelength,
            size=a.size)
    end
end

function -(a::MonoLightField2D,b::MonoLightField2D)
    if ! (a.wavelength≈b.wavelength)
        throw(DimensionMismatch("Wavelengths $(a.wavelength) and $(b.wavelength) are not equal."))
    elseif ! all(a.size .≈ b.size)
        throw(DimensionMismatch("Physical size $(a.size) and $(b.size) are not equal."))
    else
        MonoLightField2D(a.data .- b.data,
            wavelength=a.wavelength,
            size=a.size)
    end
end

*(r::Number, x::MonoLightField2D) =
    MonoLightField2D(r .* x.data,
        wavelength = x.wavelength,
        size = x.size)

*(x::MonoLightField2D, r::Number) = r*x

/(x::MonoLightField2D, r::Number) = (1/r)*x

\(r::Number, x::MonoLightField2D) = (1/r)*x

function Base.isapprox(a::MonoLightField2D,b::MonoLightField2D; kwargs...)
    isapprox(a.wavelength,b.wavelength,kwargs...) &&
      all(isapprox.(a.size,b.size,kwargs...)) &&
      isapprox(a.data,b.data,kwargs...)
end

@recipe function plot(light::MonoLightField2D)
    (m, n) = size(light.data)
    x = ((-1/2):(1/m):(1/2-1/m))*ruval(light.size[1])
    y = ((-1/2):(1/n):(1/2-1/n))*ruval(light.size[2])
    z = abs.(light.data)
    seriestype --> :heatmap
    aspect_ratio --> :equal
    xguide --> unit(light.size[1])
    yguide --> unit(light.size[2])
    size --> (m, n)
    x, y, z
end
