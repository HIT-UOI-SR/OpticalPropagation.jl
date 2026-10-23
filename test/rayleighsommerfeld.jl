# Explicit spatial quadrature: no FFT, lag indexing, or phase factorisation
# from the implementation. Compare complex fields, including absolute phase.
function rs_quadrature(U, d, λ, lx, ly; wx=ones(size(U, 1)), wy=ones(size(U, 2)),
                       output_size=(lx, ly), output_samples=size(U))
    n, m = size(U)
    dx, dy, k = lx/n, ly/m, 2*(pi/λ)
    xs = [-lx/2 + (i-1)*dx for i in 1:n]
    ys = [-ly/2 + (j-1)*dy for j in 1:m]
    nx, ny = output_samples
    lxo, lyo = output_size
    xo = [-lxo/2 + (i-1)*lxo/nx for i in 1:nx]
    yo = [-lyo/2 + (j-1)*lyo/ny for j in 1:ny]
    result = zeros(ComplexF64, nx, ny)
    for j in 1:ny, i in 1:nx, v in 1:m, u in 1:n
        r = sqrt((xo[i]-xs[u])^2 + (yo[j]-ys[v])^2 + d^2)
        result[i, j] += U[u, v]*wx[u]*wy[v]*d/(2pi) *
                        exp(im*k*r)*(1-im*k*r)/r^3 * dx*dy
    end
    result
end

@testset "Rayleigh–Sommerfeld diffraction" begin
    @testset "Full integral on $dims" for dims in
            ((1, 1), (1, 6), (5, 1), (3, 5), (4, 6), (4, 5), (5, 4), (8, 7))
        n, m = dims
        U = [complex(i-2j, i*j/3) for i in 1:n, j in 1:m]
        original = copy(U)
        # Subwavelength distance tests the usually omitted 1/r term.
        for d in (0.004, 0.7)
            expected = rs_quadrature(U, d, 0.03, 0.2, 0.15)
            @test rayleighsommerfeld(U, d, 0.03, 0.2, 0.15) ≈ expected rtol=3e-13
            @test rayleighsommerfeld(view(U, :, :), d, 0.03, 0.2, 0.15) ≈ expected rtol=3e-13
            @test rayleighsommerfeld(real.(U), d, 0.03, 0.2, 0.15) ≈
                  rs_quadrature(real.(U), d, 0.03, 0.2, 0.15) rtol=3e-13
        end
        @test U == original
    end

    @testset "Linear convolution at aperture edges" begin
        U = zeros(6, 5)
        U[1, end] = 1
        @test rayleighsommerfeld(U, 0.05, 0.04, 0.4, 0.3) ≈
              rs_quadrature(U, 0.05, 0.04, 0.4, 0.3) rtol=2e-13
        # Reversing both axes reverses the result, without wrapping an edge.
        @test rayleighsommerfeld(reverse(U), 0.05, 0.04, 0.4, 0.3) ≈
              reverse(rayleighsommerfeld(U, 0.05, 0.04, 0.4, 0.3)) rtol=2e-13
        @test rayleighsommerfeld(zeros(4, 3), 0.7, 0.03, 0.2, 0.15) == zeros(4, 3)
    end

    @testset "Simpson integration on $dims" for dims in ((3, 5), (7, 3))
        n, m = dims
        U = [complex(i-2j, i*j/3) for i in 1:n, j in 1:m]
        # Accumulate independent three-point panels, sharing their endpoints.
        wx, wy = zeros(n), zeros(m)
        for i in 1:2:n-2
            wx[i:i+2] .+= [1, 4, 1]/3
        end
        for j in 1:2:m-2
            wy[j:j+2] .+= [1, 4, 1]/3
        end
        @test rayleighsommerfeld(U, 0.13, 0.03, 0.2, 0.15; quadrature=:simpson) ≈
              rs_quadrature(U, 0.13, 0.03, 0.2, 0.15; wx=wx, wy=wy) rtol=3e-13
    end

    @testset "Analytic circular aperture on axis" begin
        # Integrating h over a disk gives exp(ikd) - d/R*exp(ikR).
        # This is a nonparaxial case with an appreciable near-field term.
        n, width, radius, λ, d = 256, 0.16, 0.06, 0.03, 0.015
        xs = [-width/2 + i*width/n for i in 0:n-1]
        U = [Float64(x^2+y^2 <= radius^2) for x in xs, y in xs]
        R = hypot(d, radius)
        expected = exp(2im*pi*d/λ) - d/R*exp(2im*pi*R/λ)
        actual = rayleighsommerfeld(U, d, λ, width, width)[n÷2+1, n÷2+1]
        # Pixelating the circular boundary dominates this quadrature error.
        @test actual ≈ expected rtol=0.003
    end

    @testset "Long axial path retains transverse phase" begin
        U = [1 2im; -3 4]
        # Float64 rounds sqrt(d²+ρ²) to d at this scale. Use high precision
        # for the reference so losing that small difference is detectable.
        expected = setprecision(256) do
            rs_quadrature(U, big"1e8", big"1e-8", big"1.0", big"0.8")
        end
        @test rayleighsommerfeld(U, 1e8, 1e-8, 1.0, 0.8) ≈ expected rtol=2e-14
    end

    @testset "Units, translated windows, functors and zero distance" begin
        U = [complex(i-j, i+2j) for i in 1:5, j in 1:3]
        field = MonoLightField2D(U; wavelength=30u"mm", size=(20u"cm", 150u"mm"),
                                center=(120u"mm", -7u"cm"))
        for quadrature in (:rectangle, :simpson)
            result = rayleighsommerfeld(field, 13u"cm"; quadrature=quadrature)
            @test parent(result) ≈ rayleighsommerfeld(U, 0.13, 0.03, 0.2, 0.15;
                                                     quadrature=quadrature)
            @test result.wavelength == field.wavelength
            @test physicalsize(result) == physicalsize(field)
            @test coordinates(result) == coordinates(field)
            @test size(result) == size(field)
            @test rayleighsommerfeld(130u"mm"; quadrature=quadrature)(field) ≈ result
            identity = rayleighsommerfeld(field, 0u"mm"; quadrature=quadrature)
            @test identity ≈ field
            @test parent(identity) !== parent(field)
            @test parent(identity) isa Matrix{ComplexF64}
        end
        identity = rayleighsommerfeld(U, 0, 0.03, 0.2, 0.15)
        @test identity == U
        @test identity !== U
        @test identity isa Matrix{ComplexF64}
    end

    @testset "Input domains" begin
        U = ones(3, 5)
        for d in (-0.1, Inf, -Inf, NaN)
            @test_throws ArgumentError rayleighsommerfeld(U, d, 0.03, 0.2, 0.15)
        end
        for invalid in (0.0, -1.0, Inf, NaN), d in (0.0, 0.13)
            @test_throws ArgumentError rayleighsommerfeld(U, d, invalid, 0.2, 0.15)
            @test_throws ArgumentError rayleighsommerfeld(U, d, 0.03, invalid, 0.15)
            @test_throws ArgumentError rayleighsommerfeld(U, d, 0.03, 0.2, invalid)
        end
        for d in (0.0, 0.13)
            @test_throws ArgumentError rayleighsommerfeld(zeros(0, 3), d, 0.03, 0.2, 0.15)
            @test_throws ArgumentError rayleighsommerfeld(zeros(3, 0), d, 0.03, 0.2, 0.15)
            @test_throws ArgumentError rayleighsommerfeld(U, d, 0.03, 0.2, 0.15; quadrature=:unknown)
            for dims in ((1, 3), (3, 1), (3, 4), (4, 3), (4, 4))
                @test_throws ArgumentError rayleighsommerfeld(ones(dims), d, 0.03, 0.2, 0.15;
                                                            quadrature=:simpson)
            end
        end
    end
end
