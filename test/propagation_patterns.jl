@testset "Canonical propagation patterns" begin
    @testset "Infinite plane wave" begin
        input = MonoLightField2D(fill(1 + 0.5im, 32, 25);
            wavelength=632.8u"nm", size=(1.2u"mm", 1u"mm"))
        λ = ustrip(u"m", input.wavelength)
        lx, ly = ustrip.(u"m", physicalsize(input))
        # A periodic window represents an unbounded plane wave without edges.
        for propagate in (angularspectrum, fresnel2), d in (-1u"cm", 1u"cm")
            phase = exp(2im*pi*Float64(d/input.wavelength))
            expected = phase * input
            @test propagate(input, d) ≈ expected rtol=5e-11
            @test propagate(parent(input), ustrip(u"m", d), λ, lx, ly) ≈
                  parent(expected) rtol=5e-11
        end
    end

    @testset "Gaussian beam" begin
        λ, waist = 500u"nm", 100u"μm"
        input = MonoLightField2D(zeros(257, 256);
            wavelength=λ, size=(2.4u"mm", 2.2u"mm"))
        xs, ys = coordinates(input)
        parent(input) .= [exp(-Float64((x^2 + y^2)/waist^2)) for x in xs, y in ys]
        rayleigh = pi*waist^2/λ
        power = sum(abs2, input)*prod(samplepitch(input))

        # The window contains the beam through two Rayleigh lengths. The
        # nonparaxial methods approach this solution as λ/waist → 0;
        # λ/waist = 0.005 allows a 2e-5 relative field error here.
        methods = ((angularspectrum, 2e-5), (fresnel2, 1e-9),
                   (fresnel1, 1e-9), (fresnelkirchhoff, 2e-5))
        @testset "$propagate, z/zR = $factor" for (propagate, tolerance) in methods,
                                                factor in (0.5, 1.0, 2.0)
            d = factor*rayleigh
            output = propagate(input, d)
            xo, yo = coordinates(output)
            q = 1 + im*Float64(d/rayleigh)
            carrier = exp(2im*pi*Float64(d/λ))
            # Complex amplitude checks spreading, wavefront curvature and Gouy phase.
            expected = [carrier/q * exp(-Float64((x^2 + y^2)/waist^2)/q)
                        for x in xo, y in yo]
            @test parent(output) ≈ expected rtol=tolerance
            @test sum(abs2, output)*prod(samplepitch(output)) ≈ power rtol=tolerance
        end
    end

    @testset "Talbot effect" begin
        λ, period = 500u"nm", 250u"μm"
        input = MonoLightField2D(zeros(96, 72);
            wavelength=λ, size=(4period, 3period))
        xs, ys = coordinates(input)
        gx(x) = 1 + 0.4*cospi(2*Float64(x/period)) + 0.2*sinpi(4*Float64(x/period))
        gy(y) = 1 + 0.3*sinpi(2*Float64(y/period))
        parent(input) .= [gx(x)*gy(y) for x in xs, y in ys]
        talbot = 2period^2/λ
        data = parent(input)
        # There are 24 samples per period on both axes; half a period is 12.
        shifted_x = circshift(data, (12, 0))
        shifted_y = circshift(data, (0, 12))
        shifted_xy = circshift(data, (12, 12))
        # At zT/4 each axis produces (1-i)/2 times the original plus
        # (1+i)/2 times its half-period shift. Expand their product in 2D.
        fractional = -0.5im*data + 0.5*(shifted_x + shifted_y) + 0.5im*shifted_xy
        images = ((1.0, data), (0.5, shifted_xy), (0.25, fractional))

        # A smooth grating has only low diffraction orders, so λ/period = 0.002
        # keeps the angular-spectrum result within 5e-5 of paraxial Talbot images.
        @testset "$propagate, z/zT = $fraction" for (propagate, tolerance) in
                ((fresnel2, 1e-9), (angularspectrum, 5e-5)), (fraction, image) in images
            d = fraction*talbot
            output = propagate(input, d)
            carrier = exp(2im*pi*Float64(d/λ))
            @test parent(output) ≈ carrier*image rtol=tolerance
        end
    end
end
