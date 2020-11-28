const LengthType=typeof(1.0u"m");


"""
    MonoLightField2D

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
struct MonoLightField2D <: AbstractArray{ComplexF64,2}
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

Base.size(a::MonoLightField2D) = size(a.data)
Base.getindex(a::MonoLightField2D, inds::Vararg{Int,2}) = a.data[inds...]
Base.setindex!(a::MonoLightField2D, val, inds::Vararg{Int,2}) = a.data[inds...] = val

Base.BroadcastStyle(::Type{MonoLightField2D}) = Broadcast.ArrayStyle{MonoLightField2D}()
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{MonoLightField2D}}, ::Type{ElType}) where ElType
    # Scan the inputs for the MonoLightField2D:
    A = find_monolightfield2d(bc)
    MonoLightField2D(A, data = similar(Array{ElType}, axes(bc)))
end

find_monolightfield2d(bc::Base.Broadcast.Broadcasted) = find_monolightfield2d(bc.args)
find_monolightfield2d(::Tuple{}) = nothing
find_monolightfield2d(args::Tuple) = find_monolightfield2d(find_monolightfield2d(args[1]), Base.tail(args))
find_monolightfield2d(x) = x
find_monolightfield2d(a::MonoLightField2D, rest) = a
find_monolightfield2d(::Any, rest) = find_monolightfield2d(rest)

@recipe function plot(light::MonoLightField2D; datafunction=abs, unit = u"mm")
    (m, n) = size(light.data)
    x = ((-1/2):(1/m):(1/2-1/m))*uconvert(Unitful.NoUnits, light.size[1]/unit)
    y = ((-1/2):(1/n):(1/2-1/n))*uconvert(Unitful.NoUnits, light.size[2]/unit)
    z = datafunction.(light.data)
    seriestype --> :heatmap
    aspect_ratio --> :equal
    xguide --> unit
    yguide --> unit
    size --> (m, n)
    x, y, z
end
