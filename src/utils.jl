function auval(x::Unitful.Length) # get absolute numerical value of quantity
    unit = u"m"
    convert(Float64, x / unit)
end

function ruval(x::Unitful.Length) # get relative numerical value of quantity
    uconvert(Unitful.NoUnits, x / unit(x))
end
