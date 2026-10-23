# Optical Propagation Methods

Here are some useful optical propagation methods.

```@docs
    angularspectrum
    fresnel1
    fresnel2
    fresnelkirchhoff
```

## Fresnel–Kirchhoff propagation

`fresnelkirchhoff` evaluates a spatial zone-mask convolution using zero-padded
FFTs. For normally incident illumination between parallel planes, the kernel
retains the spherical propagation distance and Kirchhoff obliquity factor:

```math
h(x,y) = \frac{e^{iks}}{2i\lambda s}\left(1+\frac{d}{s}\right),
\qquad s=\sqrt{d^2+x^2+y^2}, \qquad k=\frac{2\pi}{\lambda}.
```

The input is zero outside its sampled aperture, and the output covers the same
window. Sample coordinates follow the existing light-field convention:
`x[i] = -lx/2 + (i-1)*lx/n`, and likewise for `y`.
Sampling must resolve the phase of the kernel. This finite-aperture
integral requires sufficient sampling as the distance approaches zero; the
special case `d == 0` returns a copy. It does not impose periodic boundaries.

```@example kirchhoff
using OpticalPropagation, Unitful

aperture = MonoLightField2D(ones(32, 24),
    wavelength=633u"nm", size=(0.1u"mm", 0.08u"mm"))
diffracted = fresnelkirchhoff(aperture, 20u"mm")
size(diffracted), diffracted.size
```

For a spherical reference surface, set `radius` to the signed curvature radius
of the incident wave. Positive values describe divergence; negative values
describe convergence. Data contain the actual complex amplitude on the cap.
The output cap is concentric and its projected widths scale by `1 + d/radius`.
Propagation to or through the common centre is unsupported.

```@example kirchhoff
spherical = fresnelkirchhoff(aperture, 20u"mm"; radius=100u"mm")
spherical.size
```

The finite-radius mode implements the tangent-grid approximation motivated by
Hudson's Eqs. (1)–(5), with the surface-area Jacobian applied to the input.
It computes the zone mask from exact on-axis spherical geometry; translating
this mask across the projected grid approximates rotation on the reference
sphere. Use small caps relative to the radius. This grid approximation is not
an exact integral between arbitrary spherical points.

For displacement radius ``\rho=\sqrt{x^2+y^2}``, input radius ``R``, and
``c=\sqrt{1-\rho^2/R^2}``, the kernel uses

```math
z_s=R(1-c), \qquad
s=\sqrt{\rho^2+(d+z_s)^2}, \qquad
\cos\chi=\frac{dc-z_s}{s}, \qquad
h=\frac{e^{iks}}{2i\lambda s}(1+\cos\chi).
```

The quadrature weight at each input sample is ``\Delta x\Delta y/c`` using
that sample's distance from the optical axis. In the limit ``|R|\to\infty``,
these expressions reduce to the parallel-plane kernel above. The constant
reference-wave factor in Hudson's Eq. (1) is already contained in the input
amplitude; it is not multiplied into the result a second time.

Reference: Jerome A. Hudson, [“Fresnel-Kirchhoff diffraction in optical systems:
an approximate computational algorithm”](https://doi.org/10.1364/AO.23.002292),
*Applied Optics* **23**, 2292–2295 (1984). The implementation uses the derivation
in the [public publisher preview](https://www.researchgate.net/publication/5638727_Fresnel-Kirchhoff_diffraction_in_optical_systems_an_approximate_computational_algorithm).
