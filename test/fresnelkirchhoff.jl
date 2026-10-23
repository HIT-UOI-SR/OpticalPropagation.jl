# Independent quadrature in Cartesian coordinates. For a finite radius this
# integrates the actual spherical caps, without Hudson's translation approximation.
function kirchhoff_quadrature(U, d, λ, lx, ly; radius=Inf)
    n, m = size(U)
    dx, dy = lx/n, ly/m
    xs = [-lx/2 + (i - 1)*dx for i in 1:n]
    ys = [-ly/2 + (j - 1)*dy for j in 1:m]
    result = zeros(ComplexF64, n, m)
    scale = 1 + d/radius
    for j in 1:m, i in 1:n
        co = sqrt(1 - (xs[i]^2 + ys[j]^2)/radius^2)
        zo = isinf(radius) ? d : -radius + (radius + d)*co
        for v in 1:m, u in 1:n
            ci = sqrt(1 - (xs[u]^2 + ys[v]^2)/radius^2)
            zi = isinf(radius) ? 0.0 : radius*(ci - 1)
            displacement = (scale*xs[i] - xs[u], scale*ys[j] - ys[v], zo - zi)
            s = sqrt(sum(abs2, displacement))
            normal = (xs[u]/radius, ys[v]/radius, ci)
            obliquity = sum(normal[k]*displacement[k] for k in 1:3)/s
            result[i, j] += U[u, v]*exp(2im*pi*s/λ)*(1 + obliquity)*dx*dy/(2im*λ*s*ci)
        end
    end
    result
end

@testset "Fresnel–Kirchhoff diffraction" begin
    @testset "Planar integral on $dims" for dims in [(1, 1), (1, 4), (3, 1), (3, 5), (4, 6), (4, 5)]
        n, m = dims
        # Asymmetric, deterministic amplitudes expose axis swaps and FFT shifts.
        U = [complex(i - 2j, i*j/3) for i in 1:n, j in 1:m]
        original = copy(U)
        expected = kirchhoff_quadrature(U, 0.7, 0.03, 0.2, 0.15)
        @test fresnelkirchhoff(U, 0.7, 0.03, 0.2, 0.15) ≈ expected rtol=2e-13
        @test U == original
        real_input = real.(U)
        @test fresnelkirchhoff(real_input, 0.7, 0.03, 0.2, 0.15) ≈
              kirchhoff_quadrature(real_input, 0.7, 0.03, 0.2, 0.15) rtol=2e-13
        @test fresnelkirchhoff(view(U, :, :), 0.7, 0.03, 0.2, 0.15) ≈ expected rtol=2e-13
    end

    @testset "Edge impulse has no periodic wraparound" begin
        U = zeros(5, 6)
        U[1, end] = 1
        expected = kirchhoff_quadrature(U, 0.1, 0.04, 0.4, 0.3)
        @test fresnelkirchhoff(U, 0.1, 0.04, 0.4, 0.3) ≈ expected rtol=2e-13
    end

    @testset "Spherical caps, radius = $radius, size = $dims" for radius in (2.0, -2.0), dims in ((6, 4), (5, 3))
        n, m = dims
        U = [complex(i - j/2, i*j/5) for i in 1:n, j in 1:m]
        expected = kirchhoff_quadrature(U, 0.5, 0.03, 0.12, 0.08; radius=radius)
        actual = fresnelkirchhoff(U, 0.5, 0.03, 0.12, 0.08; radius=radius)
        # At the optical axis the zone-mask geometry equals the exact geometry.
        if iseven(n) && iseven(m)
            @test actual[n÷2 + 1, m÷2 + 1] ≈ expected[n÷2 + 1, m÷2 + 1] rtol=2e-12
        end
        # Away from the axis this is explicitly a small-cap approximation.
        @test actual ≈ expected rtol=0.01
    end

    @testset "Planar limit and analytic paraxial Gaussian" begin
        U = [complex(i, -j) for i in 1:3, j in 1:4]
        planar = fresnelkirchhoff(U, 0.5, 0.03, 0.12, 0.08)
        @test fresnelkirchhoff(U, 0.5, 0.03, 0.12, 0.08; radius=1e12) ≈ planar rtol=1e-10

        n, m = 129, 128
        lx, ly, λ, d, waist = 2.0, 1.8, 0.01, 10.0, 0.2
        xs = [-lx/2 + (i - 1)*lx/n for i in 1:n]
        ys = [-ly/2 + (j - 1)*ly/m for j in 1:m]
        gaussian = [exp(-(x^2 + y^2)/waist^2) for x in xs, y in ys]
        q = 1 + im*d*λ/(pi*waist^2)
        expected = [exp(2im*pi*d/λ)/q*exp(-(x^2 + y^2)/(waist^2*q)) for x in xs, y in ys]
        @test fresnelkirchhoff(gaussian, d, λ, lx, ly) ≈ expected rtol=0.002
    end

    @testset "Units, output size, functors and zero distance" begin
        U = [1 2im 3; -4 5 6im]
        field = MonoLightField2D(U, wavelength=500u"nm", size=(2u"mm", 3u"mm"))
        for radius in (Inf*u"m", 20u"cm", -20u"cm")
            result = fresnelkirchhoff(field, 10u"mm"; radius=radius)
            @test result.data ≈ fresnelkirchhoff(U, 0.01, 500e-9, 0.002, 0.003;
                                               radius=OpticalPropagation.auval(radius))
            @test result.wavelength == field.wavelength
            @test all(result.size .≈ (1 + 10u"mm"/radius) .* field.size)
            @test size(result) == size(field)
            @test fresnelkirchhoff(1u"cm"; radius=radius)(field) ≈ result
            identity = fresnelkirchhoff(field, 0u"mm"; radius=radius)
            @test identity ≈ field
            @test identity.data !== field.data
        end
        @test fresnelkirchhoff(U, 0.0, 500e-9, 0.002, 0.003) == U
        @test fresnelkirchhoff(U, 0.0, 500e-9, 0.002, 0.003) !== U
        @test fresnelkirchhoff(zeros(3, 4), 0.5, 0.03, 0.12, 0.08) == zeros(3, 4)
    end

    @testset "Invalid geometry" begin
        U = ones(3, 4)
        for d in (-1.0, Inf, NaN)
            @test_throws ArgumentError fresnelkirchhoff(U, d, 0.03, 0.12, 0.08)
        end
        for invalid in (0.0, -1.0, Inf, NaN)
            @test_throws ArgumentError fresnelkirchhoff(U, 0.5, invalid, 0.12, 0.08)
            @test_throws ArgumentError fresnelkirchhoff(U, 0.5, 0.03, invalid, 0.08)
            @test_throws ArgumentError fresnelkirchhoff(U, 0.5, 0.03, 0.12, invalid)
        end
        for radius in (0.0, -Inf, NaN, -0.5, -0.25, 0.05)
            @test_throws ArgumentError fresnelkirchhoff(U, 0.5, 0.03, 0.12, 0.08; radius=radius)
        end
        @test_throws ArgumentError fresnelkirchhoff(zeros(0, 4), 0.5, 0.03, 0.12, 0.08)
        @test_throws ArgumentError fresnelkirchhoff(zeros(3, 0), 0.5, 0.03, 0.12, 0.08)
        @test_throws ArgumentError fresnelkirchhoff(ones(1, 1), 0.5, 0.03, 0.12, 0.08; radius=0.05)
    end
end
