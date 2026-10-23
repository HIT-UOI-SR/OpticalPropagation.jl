const LengthType=typeof(1.0u"m");


"""
    MonoLightField2D

`MonoLightField2D` is used to represent a 2-dimensional monochromatic light field.

    MonoLightField2D(data; wavelength, size, center=(0u"m", 0u"m"))

`size` is the physical window width, while `Base.size` gives the sample counts.
`center` is the centre of that window, not the mean of the sample coordinates.
Dimension 1 is x and dimension 2 is y. Samples start at `center - size/2`,
with pitches `size ./ Base.size(data)`; the upper window endpoint is excluded.
The public constructors copy the amplitude data into a `Matrix{ComplexF64}`.

# Examples

```jldoctest; setup = :(using OpticalPropagation)
julia> using Unitful

julia> MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
MonoLightField2D:
    wavelength: 6.328e-7 m
    physical size: (0.001 m, 0.001 m)
    complex amplitude data:
2×2 Matrix{ComplexF64}:
 1.0+0.0im  2.0+0.0im
 3.0+0.0im  4.0+0.0im
```
"""
struct MonoLightField2D <: AbstractArray{ComplexF64,2}
    data::Array{ComplexF64,2}
    wavelength::LengthType
    size::Tuple{LengthType,LengthType}
    center::Tuple{LengthType,LengthType}

    # Internal ownership transfer for newly allocated computation results.
    function MonoLightField2D(data::Matrix{ComplexF64}, wavelength::LengthType,
                             size::Tuple{LengthType,LengthType},
                             center::Tuple{LengthType,LengthType}, ::Val{:owned})
        all(isfinite, center) || throw(ArgumentError("window centre must be finite"))
        new(data, wavelength, size, center)
    end
end

_lengthpair(x::Tuple{<:Unitful.Length,<:Unitful.Length}) =
    convert(Tuple{LengthType,LengthType}, float.(x))

function MonoLightField2D(data::AbstractMatrix; wavelength::Unitful.Length,
                         size::Tuple{<:Unitful.Length,<:Unitful.Length},
                         center::Tuple{<:Unitful.Length,<:Unitful.Length}=(0u"m", 0u"m"))
    MonoLightField2D(Matrix{ComplexF64}(data), convert(LengthType, float(wavelength)),
                     _lengthpair(size), _lengthpair(center), Val(:owned))
end

function MonoLightField2D(x::MonoLightField2D; data=x.data, wavelength=x.wavelength,
                         size=x.size, center=x.center)
    MonoLightField2D(data; wavelength=wavelength, size=size, center=center)
end

function _withdata(x::MonoLightField2D, data::Matrix{ComplexF64}; size=x.size, center=x.center)
    MonoLightField2D(data, x.wavelength, _lengthpair(size), _lengthpair(center), Val(:owned))
end
_withdata(x::MonoLightField2D, data::AbstractMatrix; kwargs...) =
    _withdata(x, Matrix{ComplexF64}(data); kwargs...)

"""
    physicalsize(light::MonoLightField2D)

Return the physical window widths `(Lx, Ly)`, also available as `light.size`.
A width is the sample count times the pitch, not the distance between endpoints.
"""
physicalsize(x::MonoLightField2D) = x.size

"""
    gridcenter(light::MonoLightField2D)

Return the physical window centre `(cx, cy)`. This is not the mean of the
sample coordinates, which lie in a half-open window.
"""
gridcenter(x::MonoLightField2D) = x.center

"""
    samplepitch(light::MonoLightField2D)

Return the sampling pitches `(Δx, Δy)` with length units. Empty dimensions
have no defined pitch and cause an `ArgumentError`.
"""
function samplepitch(x::MonoLightField2D)
    all(>(0), size(x)) || throw(ArgumentError("sample pitch is undefined for an empty dimension"))
    x.size ./ size(x)
end

"""
    coordinates(light::MonoLightField2D)
    coordinates(light::MonoLightField2D, dim)

Return unitful coordinate ranges `(x, y)`, or the range for dimension 1 or 2.
Ranges use constant storage and start at `gridcenter(light) - physicalsize(light)/2`.
`Base.axes(light)` continues to return ordinary integer indices.
"""
coordinates(x::MonoLightField2D) = (coordinates(x, 1), coordinates(x, 2))
function coordinates(x::MonoLightField2D, dim::Integer)
    dim in (1, 2) || throw(ArgumentError("a light field has only dimensions 1 and 2"))
    n = size(x, dim)
    n > 0 || throw(ArgumentError("coordinates are undefined for an empty dimension"))
    range(x.center[dim] - x.size[dim]/2; step=x.size[dim]/n, length=n)
end

"""
    crop(light::MonoLightField2D, I, J)

Copy a nonempty rectangular region into a new light field, preserving its
wavelength, sampling pitches and physical sample positions. `I` and `J` must
be contiguous integer unit ranges or `:`. Ordinary indexing retains Julia's
array behaviour and does not attach physical coordinates to slices.
"""
function crop(x::MonoLightField2D, I::Union{Colon,AbstractUnitRange{<:Integer}},
              J::Union{Colon,AbstractUnitRange{<:Integer}})
    rows = I isa Colon ? axes(x, 1) : I
    cols = J isa Colon ? axes(x, 2) : J
    checkbounds(x.data, rows, cols)
    (isempty(rows) || isempty(cols)) && throw(ArgumentError("a crop must be nonempty"))
    dx, dy = samplepitch(x)
    widths = (length(rows)*dx, length(cols)*dy)
    center = (coordinates(x, 1)[first(rows)] + widths[1]/2,
              coordinates(x, 2)[first(cols)] + widths[2]/2)
    _withdata(x, x.data[rows, cols]; size=widths, center=center)
end

function Base.show(io::IO, x::MonoLightField2D)
    print(io, "MonoLightField2D(")
    show(io, x.data)
    print(io, ", wavelength=")
    show(io, x.wavelength)
    print(io, ", size=")
    show(io, x.size)
    if any(!iszero, x.center)
        print(io, ", center=")
        show(io, x.center)
    end
    print(io, ")")
end

function Base.show(io::IO, mime::MIME"text/plain", x::MonoLightField2D)
    print(io,"MonoLightField2D:")
    print(io,"\n    wavelength: ")
    show(io, mime, x.wavelength)
    print(io, "\n    physical size: ")
    show(io, mime, x.size)
    if any(!iszero, x.center)
        print(io, "\n    window centre: ")
        show(io, mime, x.center)
    end
    print(io,"\n    complex amplitude data:\n")
    show(io, mime, x.data)
end

function _samecenter(a::MonoLightField2D, b::MonoLightField2D)
    # Compare displacement at the sampling scale, not relative to a potentially
    # very large absolute centre coordinate.
    all(1:2) do dim
        pitch = max(abs(a.size[dim]), abs(b.size[dim])) / max(size(a, dim), 1)
        isapprox(a.center[dim], b.center[dim]; rtol=0, atol=sqrt(eps(Float64))*pitch)
    end
end

function _check_compatible(a::MonoLightField2D, b::MonoLightField2D)
    if size(a) != size(b)
        throw(DimensionMismatch("Sample counts $(size(a)) and $(size(b)) are not equal."))
    elseif ! (a.wavelength≈b.wavelength)
        throw(DimensionMismatch("Wavelengths $(a.wavelength) and $(b.wavelength) are not equal."))
    elseif ! all(a.size .≈ b.size)
        throw(DimensionMismatch("Physical size $(a.size) and $(b.size) are not equal."))
    elseif !_samecenter(a, b)
        throw(DimensionMismatch("Window centres $(a.center) and $(b.center) are not equal."))
    end
    nothing
end

function +(a::MonoLightField2D,b::MonoLightField2D)
    _check_compatible(a, b)
    _withdata(a, a.data .+ b.data)
end

function -(a::MonoLightField2D,b::MonoLightField2D)
    _check_compatible(a, b)
    _withdata(a, a.data .- b.data)
end

*(r::Number, x::MonoLightField2D) =
    _withdata(x, r .* x.data)

*(x::MonoLightField2D, r::Number) = r*x

/(x::MonoLightField2D, r::Number) = (1/r)*x

\(r::Number, x::MonoLightField2D) = (1/r)*x

function Base.isapprox(a::MonoLightField2D,b::MonoLightField2D; kwargs...)
    size(a) == size(b) && _samecenter(a, b) &&
      isapprox(a.wavelength,b.wavelength; kwargs...) &&
      all(isapprox.(a.size,b.size; kwargs...)) &&
      isapprox(a.data,b.data; kwargs...)
end

Base.size(a::MonoLightField2D) = size(a.data)
Base.parent(a::MonoLightField2D) = a.data
Base.copy(a::MonoLightField2D) = _withdata(a, copy(a.data))
Base.getindex(a::MonoLightField2D, inds::Vararg{Int,2}) = a.data[inds...]
Base.setindex!(a::MonoLightField2D, val, inds::Vararg{Int,2}) = a.data[inds...] = val

Base.BroadcastStyle(::Type{MonoLightField2D}) = Broadcast.ArrayStyle{MonoLightField2D}()
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{MonoLightField2D}}, ::Type{ElType}) where ElType
    A = find_monolightfield2d(bc)
    _check_broadcast_fields(A, bc)
    axes(bc) == axes(A) || throw(DimensionMismatch("broadcasting cannot expand a physical grid"))
    _withdata(A, Matrix{ComplexF64}(undef, size(A)))
end

find_monolightfield2d(bc::Base.Broadcast.Broadcasted) = find_monolightfield2d(bc.args)
find_monolightfield2d(::Tuple{}) = nothing
find_monolightfield2d(args::Tuple) = find_monolightfield2d(find_monolightfield2d(args[1]), Base.tail(args))
find_monolightfield2d(x) = nothing
find_monolightfield2d(x::MonoLightField2D) = x
find_monolightfield2d(a::MonoLightField2D, rest) = a
find_monolightfield2d(::Any, rest) = find_monolightfield2d(rest)

_check_broadcast_fields(ref, x) = nothing
_check_broadcast_fields(ref, x::MonoLightField2D) = _check_compatible(ref, x)
_check_broadcast_fields(ref, bc::Broadcast.Broadcasted) = _check_broadcast_fields(ref, bc.args)
_check_broadcast_fields(ref, ::Tuple{}) = nothing
function _check_broadcast_fields(ref, args::Tuple)
    _check_broadcast_fields(ref, first(args))
    _check_broadcast_fields(ref, Base.tail(args))
end

# After one metadata check, let normal matrix broadcasting handle fusion and
# aliasing. No grid checks or coordinate calculations enter the pixel loop.
_matrix_broadcast_arg(x) = x
_matrix_broadcast_arg(x::MonoLightField2D) = x.data
_matrix_broadcast_arg(bc::Broadcast.Broadcasted) =
    Broadcast.Broadcasted{Nothing}(bc.f, map(_matrix_broadcast_arg, bc.args), bc.axes)

function Base.copyto!(dest::AbstractArray,
                      bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{MonoLightField2D}})
    ref = dest isa MonoLightField2D ? dest : find_monolightfield2d(bc)
    _check_broadcast_fields(ref, bc)
    axes(dest) == axes(ref) || throw(DimensionMismatch("broadcasting cannot expand a physical grid"))
    copyto!(dest isa MonoLightField2D ? dest.data : dest, _matrix_broadcast_arg(bc))
    dest
end

function Base.copyto!(dest::MonoLightField2D, bc::Broadcast.Broadcasted{Nothing})
    _check_broadcast_fields(dest, bc)
    copyto!(dest.data, _matrix_broadcast_arg(bc))
    dest
end

@recipe function plot(light::MonoLightField2D; datafunction=abs, unit = u"mm")
    (m, n) = size(light.data)
    xs, ys = coordinates(light)
    x = range(ustrip(unit, first(xs)); step=ustrip(unit, step(xs)), length=m)
    y = range(ustrip(unit, first(ys)); step=ustrip(unit, step(ys)), length=n)
    z = datafunction.(light.data)
    seriestype --> :heatmap
    aspect_ratio --> :equal
    xguide --> unit
    yguide --> unit
    size --> (m, n)
    x, y, z
end
