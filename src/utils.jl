function uval(x::Unitful.Length) # get numerical value of quantity
    unit = u"m"
    convert(Float64, x / unit)
end
