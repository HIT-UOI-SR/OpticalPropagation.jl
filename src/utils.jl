function auval(x::Unitful.Length) # get absolute numerical value of quantity
    unit = u"m"
    convert(Float64, x / unit)
end
