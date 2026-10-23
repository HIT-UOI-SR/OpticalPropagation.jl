# Independent references use explicit sums, not FFTs, frequency helpers or the
# chirp factorisation in the implementation. Compare complex amplitudes so that
# origin phases, transform normalisation and axis orientation are all observable.
function transfer_dft_reference(input, d, λ, lx, ly; angular=false)
    n, m = size(input)
    spectrum = [sum(input[i+1, j+1] * exp(-2im*pi*(k*i/n + l*j/m))
                    for i in 0:n-1, j in 0:m-1) for k in 0:n-1, l in 0:m-1]
    for l in 0:m-1, k in 0:n-1
        fx = (k <= (n-1)÷2 ? k : k-n)/lx
        fy = (l <= (m-1)÷2 ? l : l-m)/ly
        if angular
            longitudinal = 1 - (λ*fx)^2 - (λ*fy)^2
            transfer = longitudinal >= 0 ? exp(2im*pi*d/λ*sqrt(longitudinal)) : 0im
        else
            transfer = exp(2im*pi*d/λ) * exp(-im*pi*λ*d*(fx^2 + fy^2))
        end
        spectrum[k+1, l+1] *= transfer
    end
    [sum(spectrum[k+1, l+1] * exp(2im*pi*(k*i/n + l*j/m))
         for k in 0:n-1, l in 0:m-1)/(n*m) for i in 0:n-1, j in 0:m-1]
end

function fresnel_quadrature(input, d, λ, xs, ys, xo, yo, dx, dy)
    [exp(2im*pi*d/λ)/(im*λ*d) * dx*dy *
     sum(input[i, j] * exp(im*pi/(λ*d)*((x-xs[i])^2 + (y-ys[j])^2))
         for i in axes(input, 1), j in axes(input, 2)) for x in xo, y in yo]
end

@testset "Propagation numerics" begin
    shapes = ((1, 1), (1, 5), (5, 1), (4, 6), (5, 3), (4, 5), (5, 4))

    @testset "FFT-order plane waves: $dims" for dims in shapes
        n, m = dims
        λ, lx, ly = 0.01, 0.2, 0.15
        modes = unique([(0, 0), ((n-1)÷2, 0), (0, (m-1)÷2),
                        (-n÷2, 0), (0, -m÷2), (-n÷2, -m÷2)])
        for (kx, ky) in modes, d in (0.7, -0.7)
            input = [exp(2im*pi*(kx*i/n + ky*j/m)) for i in 0:n-1, j in 0:m-1]
            fx, fy = kx/lx, ky/ly
            expected_as = input .* exp(2im*pi*d/λ*sqrt(1-(λ*fx)^2-(λ*fy)^2))
            expected_fresnel = input .* (exp(2im*pi*d/λ)*exp(-im*pi*λ*d*(fx^2+fy^2)))
            @test angularspectrum(input, d, λ, lx, ly) ≈ expected_as rtol=3e-12
            @test fresnel2(input, d, λ, lx, ly) ≈ expected_fresnel rtol=3e-12
        end
    end

    @testset "Direct DFT on asymmetric grids: $dims" for dims in shapes
        n, m = dims
        input = [complex(i-2j, i*j/3) for i in 1:n, j in 1:m]
        original = copy(input)
        for d in (0.7, -0.7), λ in (0.01, 0.2)
            @test angularspectrum(input, d, λ, 0.2, 0.15) ≈
                transfer_dft_reference(input, d, λ, 0.2, 0.15; angular=true) rtol=2e-12
            @test fresnel2(input, d, λ, 0.2, 0.15) ≈
                transfer_dft_reference(input, d, λ, 0.2, 0.15) rtol=2e-12
        end
        @test input == original
        for propagate in (angularspectrum, fresnel2)
            @test propagate(view(input, :, :), 0.7, 0.01, 0.2, 0.15) ≈
                propagate(input, 0.7, 0.01, 0.2, 0.15)
            forward = propagate(input, 0.7, 0.01, 0.2, 0.15)
            @test propagate(forward, -0.7, 0.01, 0.2, 0.15) ≈ input rtol=2e-12
        end
    end

    @testset "Evanescent filtering and grazing modes" begin
        input = [isodd(i) ? -1.0 : 1.0 for i in 0:3, j in 0:2]
        for d in (0.7, -0.7)
            @test maximum(abs, angularspectrum(input, d, 1.0, 1.0, 1.0)) < 1e-12
            # At the cutoff the longitudinal wave number and propagation phase
            # are zero, so the grazing mode is retained without a phase change.
            @test angularspectrum(input, d, 0.5, 1.0, 1.0) ≈ input rtol=2e-12
        end
    end

    @testset "Single-FFT direct Fresnel quadrature: $dims" for dims in shapes
        n, m = dims
        λ, lx, ly = 0.03, 0.2, 0.15
        dx, dy = lx/n, ly/m
        xs = [-lx/2 + i*dx for i in 0:n-1]
        ys = [-ly/2 + j*dy for j in 0:m-1]
        input = [complex(i-2j, i*j/3) for i in 1:n, j in 1:m]
        original = copy(input)
        for d in (0.7, -0.7)
            actual, (lxo, lyo) = fresnel1(input, d, λ, lx, ly)
            @test lxo ≈ λ*abs(d)/dx
            @test lyo ≈ λ*abs(d)/dy
            @test lxo > 0 && lyo > 0
            xo = [-lxo/2 + i*lxo/n for i in 0:n-1]
            yo = [-lyo/2 + j*lyo/m for j in 0:m-1]
            expected = fresnel_quadrature(input, d, λ, xs, ys, xo, yo, dx, dy)
            @test actual ≈ expected rtol=2e-12
            @test input == original
            @test sum(abs2, actual)*(lxo/n)*(lyo/m) ≈ sum(abs2, input)*dx*dy rtol=2e-12
            reconstructed, original_widths = fresnel1(actual, -d, λ, lxo, lyo)
            @test reconstructed ≈ input rtol=2e-12
            @test all(original_widths .≈ (lx, ly))
            @test fresnel1(view(input, :, :), d, λ, lx, ly)[1] ≈ expected rtol=2e-12
            real_input = real.(input)
            @test fresnel1(real_input, d, λ, lx, ly)[1] ≈
                fresnel_quadrature(real_input, d, λ, xs, ys, xo, yo, dx, dy) rtol=2e-12
        end
    end

    @testset "Unitful coordinates and translated windows" begin
        input = MonoLightField2D([complex(i-j, i+2j) for i in 1:5, j in 1:4];
            wavelength=30u"mm", size=(20u"cm", 150u"mm"), center=(120u"mm", -7u"cm"))
        xs, ys = coordinates(input)
        dx, dy = samplepitch(input)
        for d in (70u"cm", -700u"mm")
            output = fresnel1(input, d)
            xo, yo = coordinates(output)
            expected = fresnel_quadrature(parent(input), ustrip(u"m", d), 0.03,
                ustrip.(u"m", xs), ustrip.(u"m", ys),
                ustrip.(u"m", xo), ustrip.(u"m", yo), ustrip(u"m", dx), ustrip(u"m", dy))
            @test parent(output) ≈ expected rtol=2e-12
            @test gridcenter(output) == gridcenter(input)
            @test all(>(0u"m"), physicalsize(output))
            @test all(>(0u"m"), samplepitch(output))
            @test fresnel1(-d)(output) ≈ input rtol=2e-12
        end
    end

    @testset "Analytic Gaussian, including output phase" begin
        n, m = 129, 128
        lx, ly, λ, d, waist = 2.0, 1.8, 0.01, 10.0, 0.2
        xs = [-lx/2 + i*lx/n for i in 0:n-1]
        ys = [-ly/2 + j*ly/m for j in 0:m-1]
        input = [exp(-(x^2+y^2)/waist^2) for x in xs, y in ys]
        for distance in (d, -d)
            actual, (lxo, lyo) = fresnel1(input, distance, λ, lx, ly)
            xo = [-lxo/2 + i*lxo/n for i in 0:n-1]
            yo = [-lyo/2 + j*lyo/m for j in 0:m-1]
            q = 1 + im*distance*λ/(pi*waist^2)
            expected = [exp(2im*pi*distance/λ)/q * exp(-(x^2+y^2)/(waist^2*q)) for x in xo, y in yo]
            @test actual ≈ expected rtol=1e-7
        end
    end

    @testset "Zero distance and input domains" begin
        # The zero-distance identity also retains components that the angular
        # spectrum implementation filters for any nonzero distance.
        input = [1 2im 3; -4 5 6im]
        for propagate in (angularspectrum, fresnel2, fresnel1)
            result = propagate(input, 0.0, 1.0, 0.2, 0.15)
            data = propagate === fresnel1 ? result[1] : result
            @test data == input
            @test data !== input
            @test data isa Matrix{ComplexF64}
            if propagate === fresnel1
                @test result[2] == (0.2, 0.15)
            end
            field = MonoLightField2D(input; wavelength=633u"nm", size=(2u"mm", 3u"mm"),
                                     center=(5u"mm", -6u"mm"))
            identity = propagate(field, 0u"mm")
            @test identity ≈ field
            @test parent(identity) !== parent(field)
            @test gridcenter(identity) == gridcenter(field)
            for invalid in (Inf, -Inf, NaN)
                @test_throws ArgumentError propagate(input, invalid, 0.03, 0.2, 0.15)
            end
            for invalid in (0.0, -1.0, Inf, NaN)
                @test_throws ArgumentError propagate(input, 0.7, invalid, 0.2, 0.15)
                @test_throws ArgumentError propagate(input, 0.7, 0.03, invalid, 0.15)
                @test_throws ArgumentError propagate(input, 0.7, 0.03, 0.2, invalid)
            end
            @test_throws ArgumentError propagate(zeros(0, 3), 0.7, 0.03, 0.2, 0.15)
            @test_throws ArgumentError propagate(zeros(3, 0), 0.7, 0.03, 0.2, 0.15)
        end
    end
end
