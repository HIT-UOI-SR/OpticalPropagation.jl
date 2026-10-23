# Optical Propagation Methods

Here are some useful optical propagation methods.

```@docs
    angularspectrum
    fresnel1
    fresnel2
    fresnelkirchhoff
    rayleighsommerfeld
```

## Physical windows

The light-field interfaces preserve `gridcenter(field)` for angular-spectrum,
double-FFT Fresnel, single-FFT Fresnel and planar Fresnel–Kirchhoff propagation.
The output window follows the input centre. Single-FFT Fresnel updates the
physical widths using its sampling relation; the other planar methods preserve
the widths. Single-FFT widths are positive for either propagation direction,
so negative distances no longer reverse the output coordinate axes.

All three planar FFT methods accept finite signed distances, finite positive
wavelengths and input widths, and nonempty one-based matrices. At zero distance
they return an independent complex copy on the original grid. The matrix
interfaces retain their original signatures and return conventions.

Finite-radius Fresnel–Kirchhoff propagation currently requires a zero-centred
window. It rejects nonzero centres because the existing surface weights and
geometry checks assume coordinates about the optical axis. Off-axis spherical
propagation needs a separate numerical extension. Centred spherical windows
continue to scale by `1 + d/radius`.

## FFT frequency grids

Angular-spectrum and double-FFT Fresnel propagation apply transfer functions in
the native FFT frequency order. For an axis with `N` samples over a width `L`,
this is `fftfreq(N, N/L)`, including the zero-frequency sample for odd sizes.
The spatial coordinates still follow the existing half-open window convention.
Frequency bins and spatial sample positions have separate roles; an odd spatial
window does not imply a half-bin displacement of the frequency grid.

These two propagators impose periodic boundaries on the input window. The angular
spectrum transfer suppresses evanescent components for any nonzero distance,
including backward propagation, rather than exponentially amplifying them.
Consequently its forward/backward pair recovers only the propagating spectrum.
The exact zero-distance identity retains all input components. Components exactly
at the propagating cutoff are retained with zero longitudinal phase.

See [AbstractFFTs' frequency and transform definitions](https://juliamath.github.io/AbstractFFTs.jl/stable/api/).

## Single-FFT Fresnel propagation

For nonzero distance, `fresnel1` evaluates the discrete Fresnel integral

```math
U_o(x_o,y_o) = \frac{e^{2\pi i d/\lambda}}{i\lambda d}\Delta x\Delta y
\sum_{i,j} U_{i,j}\exp\left[\frac{i\pi}{\lambda d}
\left((x_o-x_i)^2+(y_o-y_j)^2\right)\right].
```

The output has the same sample counts, with

```math
\Delta x_o = \frac{\lambda |d|}{N_x\Delta x},\qquad
\Delta y_o = \frac{\lambda |d|}{N_y\Delta y},\qquad
L_{x,o}=N_x\Delta x_o,\quad L_{y,o}=N_y\Delta y_o.
```

Input and output coordinates both start at their window centre minus half the
corresponding width. The implementation includes the input and output quadratic
phases and the phase associated with those coordinate origins. This applies to
odd sizes and singleton dimensions as well as even grids. Both windows can be
translated by the same centre without changing the result in their local coordinates.

Positive distance uses a forward FFT; negative distance uses an unnormalised
backward FFT with the corresponding phase and integral prefactor. Output widths
remain positive in either case. The carrier phase is retained, so results can be
compared as complex amplitudes, without adjusting an arbitrary constant phase.
At zero distance the FFT sampling relation is bypassed and the input grid is kept.

The method preserves the discrete power `sum(abs2, U) * Δx * Δy` under its natural
output sampling. Physical accuracy still requires the input window and spacing
to resolve the field and quadratic phase; no automatic resampling or aliasing
correction is performed.

```@example singlefft
using OpticalPropagation, Unitful

field = MonoLightField2D(ones(5, 4),
    wavelength=633u"nm", size=(0.1u"mm", 0.08u"mm"),
    center=(1u"mm", -2u"mm"))
forward = fresnel1(field, 2u"mm")
backward = fresnel1(forward, -2u"mm")
all(>(0u"m"), physicalsize(forward)), gridcenter(forward), isapprox(backward, field)
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

## Rayleigh–Sommerfeld propagation

`rayleighsommerfeld` evaluates the first Rayleigh–Sommerfeld integral using
Shen and Wang's FFT direct integration (FFT-DI). It convolves the finite input
aperture with the full impulse response

```math
h(x,y;d) = \frac{d\,e^{ikr}}{2\pi r^2}\left(\frac{1}{r}-ik\right),
\qquad r=\sqrt{x^2+y^2+d^2},\qquad k=\frac{2\pi}{\lambda}.
```

The spherical phase, obliquity factor and near-field term are retained, without
the Fresnel approximation. For equal grids, padding to at least
`(2Nx-1, 2Ny-1)` prevents periodic
wraparound in the retained output window; FFT lengths are rounded up to products
of 2, 3 and 5 for efficiency. By default the output preserves the sample counts
and physical widths. The wavelength and window centre are always preserved;
both windows can be translated together. Distances must be finite and
nonnegative. Zero distance returns a complex copy on the original grid;
requesting a different output grid at zero distance raises an `ArgumentError`.

The default `quadrature=:rectangle` evaluates the Riemann sum with area
`Δx*Δy` for every sample. It supports rectangular arrays, odd and even counts,
and singleton dimensions. The optional `quadrature=:simpson` applies the tensor
product of `[1,4,2,…,4,1]/3` before convolution. Each dimension must then contain
an odd number of at least three samples. Simpson's rule integrates between the
first and last sample coordinates, a span of `(N-1)*Δ`; it does not add a sample
at the excluded upper endpoint of the physical window of width `N*Δ`.

```@example rayleighsommerfeld
using OpticalPropagation, Unitful

aperture = MonoLightField2D(ones(65, 49),
    wavelength=633u"nm", size=(0.1u"mm", 0.08u"mm"))
diffracted = rayleighsommerfeld(aperture, 20u"mm")
simpson = rayleighsommerfeld(20u"mm"; quadrature=:simpson)(aperture)
size(diffracted), samplepitch(simpson), gridcenter(simpson)
```

Resolve both the input field and the oscillations of the kernel, especially for
small distances or large transverse displacements. Zero padding prevents circular
convolution but does not repair undersampling or an insufficient input aperture.
Refine the sample pitch and check convergence when choosing a grid. Finite output
windows may omit diffracted power, so a propagation pair is not generally an
inverse and need not conserve the power measured inside that window.

### Independent input and output sampling

Set `output_size` and `output_samples` to select the output window widths and
sample counts independently, including different scale factors along x and y.
Each option defaults to the input value. Widths have length units for a
`MonoLightField2D` and are in metres for the matrix interface. The matrix
interface returns output amplitudes. The light-field result also stores the
selected widths.
Coordinates follow the same half-open convention on both grids:
`x[i] = center_x - Lx/2 + (i-1)*Lx/Nx`. Odd counts retain their half-sample origin
offset. The output window follows the input centre.

```@example rayleighsommerfeld
zoomed = rayleighsommerfeld(aperture, 20u"mm";
    output_size=(0.06u"mm", 0.12u"mm"), output_samples=(81, 72))
physicalsize(zoomed), size(zoomed), samplepitch(zoomed)
```

A changed grid selects Nascov and Logofătu's scaled convolution, Eqs. (14)–(22),
generalized to unequal counts and rectangular arrays. It retains this package's
`exp(+ikr)` propagation convention and the full near-field term. For each axis,
the kernel pitch and two scale factors are

```math
\delta=\max(\Delta_i,\Delta_o),\qquad
\alpha_i=\Delta_i/\delta,\qquad \alpha_o=\Delta_o/\delta.
```

The kernel is sampled on a centred, even grid with periods `Mx` and `My` that
enclose every source-to-target displacement. The input's scaled Fourier transform
is multiplied by the ordinary transform of the kernel, then transformed back at
the output scale. With zero-based input indices ``i,j`` and centred frequency
indices ``p,q``, this gives

```math
\widehat U_{p,q}=\sum_{i,j} w_{i,j} U_{i,j}
\exp\left[-2\pi i\left(
\frac{\alpha_{i,x}(i-N_{i,x}/2)p}{M_x}+
\frac{\alpha_{i,y}(j-N_{i,y}/2)q}{M_y}\right)\right],
```

```math
U_{o,a,b}=\frac{\Delta_{i,x}\Delta_{i,y}}{M_xM_y}
\sum_{p,q} H_{p,q}\widehat U_{p,q}
\exp\left[2\pi i\left(
\frac{\alpha_{o,x}(a-N_{o,x}/2)p}{M_x}+
\frac{\alpha_{o,y}(b-N_{o,y}/2)q}{M_y}\right)\right].
```

Here ``w`` contains the rectangle or Simpson input weights, and ``H`` is the DFT
of the sampled kernel ``h(s\delta_x,t\delta_y)``. The area factor uses the input
pitches. The scaled transforms use the chirp convolution from Appendix A, with
zero-padded FFTs; no dense transform matrices or direct source-to-target sums are
constructed. With comparable counts on both axes, the work scales as
``O(M_xM_y\log(M_xM_y))``. Simpson restrictions apply to the input counts only;
the output counts may be even or odd.

Unlike the equal-grid FFT-DI path, scaled convolution evaluates a periodic Fourier
interpolant of the kernel at fractional displacements. It therefore has an
interpolation error in addition to input quadrature error. Resolve the kernel on
**both** grids: a fine input grid alone does not compensate for a coarse output
grid because the larger pitch determines ``\delta``. For each axis, the phase
sampling condition from Section 4 is

```math
\frac{\delta X}{\lambda\sqrt{d^2+X^2}}<\frac12,\qquad X=M\delta/2.
```

Also resolve the near-field amplitude variations on the scale of `d`. These are
sampling guidelines; verify convergence by refining input and output pitches.
Increasing the positive integer `padding` multiplies the kernel periods, moving
the periodic boundary away from the needed displacements. It increases memory
and work and does not refine the kernel pitch; it cannot correct undersampling.
For example, add `padding=2` to compare with twice the default period on each axis.
For an unchanged output grid, the original Shen–Wang FFT-DI path is used and
`padding` has no effect after validation.

References:

1. F. Shen and A. Wang, [“Fast-Fourier-transform based numerical integration method
   for the Rayleigh-Sommerfeld diffraction formula”](https://doi.org/10.1364/AO.45.001102),
   *Applied Optics* **45**, 1102–1110 (2006).
   [Open full text](https://vtechworks.lib.vt.edu/bitstream/handle/10919/46900/ao-45-6-1102.pdf).
2. V. Nascov and P. C. Logofătu, [“Fast computation algorithm for the
   Rayleigh-Sommerfeld diffraction formula using a type of scaled convolution”](https://doi.org/10.1364/AO.48.004310),
   *Applied Optics* **48**, 4310–4319 (2009).
   [Author-provided full text](https://www.researchgate.net/publication/26711403_Fast_computation_algorithm_for_the_Rayleigh-Sommerfeld_diffraction_formula_using_a_type_of_scaled_convolution).
