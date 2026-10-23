# Independent sums over the centred Fourier basis. This deliberately uses no
# FFTs, chirps, implementation helpers, or automatic kernel-period selection.
function rs_fourier_reference(U, d, λ, Li, Lo, no, periods;
                              wx=ones(size(U, 1)), wy=ones(size(U, 2)))
    ni = size(U)
    di, do_ = Li ./ ni, Lo ./ no
    δ = max.(di, do_)
    nx, ny = periods
    ks, ls = collect(-nx÷2:nx÷2-1), collect(-ny÷2:ny÷2-1)
    xs = [-Li[1]/2 + i*di[1] for i in 0:ni[1]-1]
    ys = [-Li[2]/2 + j*di[2] for j in 0:ni[2]-1]
    xo = [-Lo[1]/2 + i*do_[1] for i in 0:no[1]-1]
    yo = [-Lo[2]/2 + j*do_[2] for j in 0:no[2]-1]
    h = [begin
        r = sqrt((s*δ[1])^2 + (t*δ[2])^2 + d^2)
        d/(2pi)*exp(2im*pi*r/λ)*(1-2im*pi*r/λ)/r^3
    end for s in ks, t in ls]
    H = [sum(h[i,j]*exp(-2im*pi*(k*ks[i]/nx + l*ls[j]/ny))
             for i in eachindex(ks), j in eachindex(ls)) for k in ks, l in ls]
    F = [sum(U[i,j]*wx[i]*wy[j]*exp(-2im*pi*(k*xs[i]/(nx*δ[1]) + l*ys[j]/(ny*δ[2])))
             for i in 1:ni[1], j in 1:ni[2]) for k in ks, l in ls]
    [prod(di)/prod(periods)*sum(F[i,j]*H[i,j]*
        exp(2im*pi*(ks[i]*x/(nx*δ[1]) + ls[j]*y/(ny*δ[2])))
        for i in eachindex(ks), j in eachindex(ls)) for x in xo, y in yo]
end

@testset "Rayleigh–Sommerfeld scaled convolution" begin
    @testset "Rectangular scaled DFT, $dims" for dims in
            ((1, 1), (1, 4), (5, 1), (3, 4), (4, 3))
        A = [complex(i-j, i+2j) for i in 1:dims[1], j in 1:dims[2]]
        for dim in 1:2, nout in (1, 6), β in (-0.073, 0.11)
            nin = size(A, dim)
            a, b = -nin/2, -nout/2
            expected = if dim == 1
                [sum(A[i+1,k]*exp(2im*pi*β*(i+a)*(j+b)) for i in 0:nin-1)
                 for j in 0:nout-1, k in axes(A,2)]
            else
                [sum(A[k,i+1]*exp(2im*pi*β*(i+a)*(j+b)) for i in 0:nin-1)
                 for k in axes(A,1), j in 0:nout-1]
            end
            @test OpticalPropagation._rs_scaled_dft(A, dim, nout, β, a, b) ≈ expected rtol=2e-13
        end
    end

    cases = (((5,4), (4,3), (0.02,0.016), (0.027,0.008), (10,8)),
             ((3,5), (4,6), (0.021,0.017), (0.012,0.025), (8,16)))
    @testset "Centred Fourier reference, $ni to $no" for (ni, no, Li, Lo, periods) in cases
        U = [complex(i-j/2, i*j/5) for i in 1:ni[1], j in 1:ni[2]]
        original = copy(U)
        for padding in (1, 2)
            expected = rs_fourier_reference(U, 0.11, 0.03, Li, Lo, no, padding .* periods)
            actual = rayleighsommerfeld(U, 0.11, 0.03, Li...;
                output_size=Lo, output_samples=no, padding=padding)
            @test size(actual) == no
            @test actual isa Matrix{ComplexF64}
            @test actual ≈ expected rtol=2e-12
            @test rayleighsommerfeld(view(U,:,:), 0.11, 0.03, Li...;
                output_size=Lo, output_samples=no, padding=padding) ≈ expected rtol=2e-12
        end
        @test U == original
        @test rayleighsommerfeld(real.(U), 0.11, 0.03, Li...; output_size=Lo, output_samples=no) ≈
              rs_fourier_reference(real.(U), 0.11, 0.03, Li, Lo, no, periods) rtol=2e-12
        if ni == (3,5)
            wx, wy = [1,4,1]/3, [1,4,2,4,1]/3
            @test rayleighsommerfeld(U, 0.11, 0.03, Li...;
                output_size=Lo, output_samples=no, quadrature=:simpson) ≈
                rs_fourier_reference(U, 0.11, 0.03, Li, Lo, no, periods; wx=wx, wy=wy) rtol=2e-12
        end
    end

    @testset "Equal pitches with different counts: $ni to $no" for (ni,no) in
            (((4,6),(8,4)), ((3,5),(7,3)), ((1,1),(3,5)), ((5,1),(3,1)))
        di = (0.01, 0.007)
        Li, Lo = ni .* di, no .* di
        U = zeros(ComplexF64, ni)
        U[1,end] = 1 + 0.5im
        # Same-parity counts give integer kernel displacements: Fourier
        # interpolation then recovers the sampled kernel exactly, even near-field.
        @test rayleighsommerfeld(U, 0.004, 0.03, Li...; output_size=Lo, output_samples=no) ≈
              rs_quadrature(U, 0.004, 0.03, Li...; output_size=Lo, output_samples=no) rtol=2e-12
    end

    @testset "Convergence to direct integration" begin
        Li, Lo = (0.02,0.016), (0.03,0.01)
        errors = Float64[]
        for (ni,no) in (((17,16),(20,13)), ((33,32),(40,26)), ((65,64),(80,52)))
            xs = [-Li[1]/2 + i*Li[1]/ni[1] for i in 0:ni[1]-1]
            ys = [-Li[2]/2 + j*Li[2]/ni[2] for j in 0:ni[2]-1]
            U = [exp(-(x^2+y^2)/0.008^2)*(1+0.3im*x/0.01-0.2im*y/0.01) for x in xs, y in ys]
            expected = rs_quadrature(U, 0.1, 0.03, Li...; output_size=Lo, output_samples=no)
            actual = rayleighsommerfeld(U, 0.1, 0.03, Li...; output_size=Lo, output_samples=no, padding=2)
            push!(errors, sqrt(sum(abs2,actual-expected)/sum(abs2,expected)))
        end
        @test errors[2] < 0.3errors[1]
        @test errors[3] < 0.3errors[2]
        @test errors[3] < 5e-7
    end

    @testset "Analytic Gaussian on an independently sampled output" begin
        ni, no = (129,128), (161,142)
        Li, Lo = (0.0012,0.0011), (0.0006,0.0009)
        waist, λ = 0.0001, 5e-7
        d = pi*waist^2/λ
        xs = [-Li[1]/2 + i*Li[1]/ni[1] for i in 0:ni[1]-1]
        ys = [-Li[2]/2 + j*Li[2]/ni[2] for j in 0:ni[2]-1]
        U = [exp(-(x^2+y^2)/waist^2) for x in xs, y in ys]
        xo = [-Lo[1]/2 + i*Lo[1]/no[1] for i in 0:no[1]-1]
        yo = [-Lo[2]/2 + j*Lo[2]/no[2] for j in 0:no[2]-1]
        q = 1 + im*d*λ/(pi*waist^2)
        expected = [cispi(2*d/λ)/q*exp(-(x^2+y^2)/(waist^2*q)) for x in xo, y in yo]
        actual = rayleighsommerfeld(U,d,λ,Li...; output_size=Lo,output_samples=no)
        # The paraxial Gaussian is the small-λ/waist limit of the full RS integral.
        @test actual ≈ expected rtol=3e-6
    end

    @testset "Output metadata, defaults and zero distance" begin
        U = [complex(i-j,i+j) for i in 1:3,j in 1:5]
        field = MonoLightField2D(U; wavelength=30u"mm", size=(21u"mm",17u"mm"),
                                center=(3u"cm",-2u"mm"))
        Lo, no = (12u"mm",25u"mm"), (4,6)
        output = rayleighsommerfeld(field, 11u"cm"; output_size=Lo, output_samples=no, padding=2)
        @test size(output) == no
        @test all(physicalsize(output) .≈ Lo)
        @test all(samplepitch(output) .≈ Lo ./ no)
        @test gridcenter(output) == gridcenter(field)
        @test output.wavelength == field.wavelength
        @test first(coordinates(output,1)) ≈ gridcenter(field)[1]-Lo[1]/2
        @test first(coordinates(output,2)) ≈ gridcenter(field)[2]-Lo[2]/2
        @test parent(output) ≈ rayleighsommerfeld(U, 0.11, 0.03, 0.021, 0.017;
            output_size=(0.012,0.025), output_samples=no, padding=2)
        @test rayleighsommerfeld(110u"mm"; output_size=Lo, output_samples=no, padding=2)(field) ≈ output
        @test size(rayleighsommerfeld(field, 11u"cm"; output_size=Lo)) == size(field)
        counts_only = rayleighsommerfeld(field, 11u"cm"; output_samples=no)
        @test physicalsize(counts_only) == physicalsize(field)
        @test size(counts_only) == no
        @test parent(rayleighsommerfeld(field, 11u"cm";
            output_size=physicalsize(field), output_samples=size(field))) ==
            parent(rayleighsommerfeld(field, 11u"cm"))
        identity = rayleighsommerfeld(field, 0u"m";
            output_size=physicalsize(field), output_samples=size(field), padding=2)
        @test identity ≈ field
        @test parent(identity) !== parent(field)
        @test rayleighsommerfeld(zeros(3,5), 0.1, 0.03, 0.02, 0.016;
            output_size=(0.03,0.01), output_samples=(4,6)) == zeros(4,6)
    end

    @testset "Invalid output grids" begin
        U = ones(3,5)
        for d in (0.0,0.1), bad in (0.0,-1.0,Inf,NaN)
            @test_throws ArgumentError rayleighsommerfeld(U,d,0.03,0.02,0.016; output_size=(bad,0.02))
            @test_throws ArgumentError rayleighsommerfeld(U,d,0.03,0.02,0.016; output_size=(0.02,bad))
        end
        for no in ((0,5),(3,0),(-1,5),(3,-1))
            @test_throws ArgumentError rayleighsommerfeld(U,0.1,0.03,0.02,0.016; output_samples=no)
        end
        for pad in (0,-1)
            @test_throws ArgumentError rayleighsommerfeld(U,0.1,0.03,0.02,0.016; padding=pad)
        end
        @test_throws ArgumentError rayleighsommerfeld(U,0.0,0.03,0.02,0.016; output_size=(0.03,0.01))
        @test_throws ArgumentError rayleighsommerfeld(U,0.0,0.03,0.02,0.016; output_samples=(4,6))
        @test_throws ArgumentError rayleighsommerfeld(ones(4,5),0.1,0.03,0.02,0.016;
            output_samples=(3,3), quadrature=:simpson)
    end
end
