### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ cd6e8890-3632-11eb-0528-e1c9df81dc75
using Plots

# ╔═╡ e85c9b60-3632-11eb-22ba-9b4355ebac76
using PlutoUI

# ╔═╡ ff6eac80-3632-11eb-0cd1-5965f4ebeb11
using Unitful

# ╔═╡ f2500fd2-3632-11eb-14f1-8d0d9d355f8e
using OpticalPropagation

# ╔═╡ 154668e0-3633-11eb-169e-bd2dbc3c558d
gauss = MonoLightField2D(
    [exp(-(x^2+y^2)/3200) for x in -200:199, y in -200:199],
    wavelength=632.8u"nm",
    size=(1u"mm", 1u"mm")
);

# ╔═╡ a345c4f0-3634-11eb-17e3-f9096b32f28e
@bind d Slider(-20:20, default=0, show_value=true)

# ╔═╡ 34d8d670-3633-11eb-3f29-15e1181b491a
plot(plot(gauss), plot(fresnel2(gauss, d*u"cm")), size=(800,400))

# ╔═╡ Cell order:
# ╠═cd6e8890-3632-11eb-0528-e1c9df81dc75
# ╠═e85c9b60-3632-11eb-22ba-9b4355ebac76
# ╠═ff6eac80-3632-11eb-0cd1-5965f4ebeb11
# ╠═f2500fd2-3632-11eb-14f1-8d0d9d355f8e
# ╠═154668e0-3633-11eb-169e-bd2dbc3c558d
# ╠═a345c4f0-3634-11eb-17e3-f9096b32f28e
# ╠═34d8d670-3633-11eb-3f29-15e1181b491a
