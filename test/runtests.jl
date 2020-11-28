using OpticalPropagation
using Test
using Unitful
using RecipesBase

@testset "OpticalPropagation.jl" begin
    @testset "MonoLightField2D Constructors" begin
        o=MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
        @test begin
            ro=MonoLightField2D(o,data=[3 1; 2 4])
            r=MonoLightField2D([3 1; 2 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            ro≈r
        end
        @test begin
            ro=MonoLightField2D(o,wavelength=1.15u"μm")
            r=MonoLightField2D([1 2; 3 4],wavelength=1.15u"μm",size=(1u"mm",1u"mm"))
            ro≈r
        end
        @test begin
            ro=MonoLightField2D(o,size=(1.5u"mm",2.0u"mm"))
            r=MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1.5u"mm",2.0u"mm"))
            ro≈r
        end
        @test begin
            ro=MonoLightField2D(o,data=[3 1; 2 4],wavelength=1.15u"μm")
            r=MonoLightField2D([3 1; 2 4],wavelength=1.15u"μm",size=(1u"mm",1u"mm"))
            ro≈r
        end
        @test begin
            ro=MonoLightField2D(o,wavelength=1.15u"μm",size=(1.5u"mm",2.0u"mm"))
            r=MonoLightField2D([1 2; 3 4],wavelength=1.15u"μm",size=(1.5u"mm",2.0u"mm"))
            ro≈r
        end
    end
    @testset "Additive Operator $op" for op in [(+), (-)]
        @test begin
            ad = [1 2; 3 4]
            bd = [1.0 -2; 3im -5]
            a=MonoLightField2D(ad,wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            b=MonoLightField2D(bd,wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            c=MonoLightField2D(op(ad, bd),wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            op(a,b)≈c
        end
        @test_throws DimensionMismatch begin
            a=MonoLightField2D([1 2; 3 4],wavelength=589.6u"nm",size=(1u"mm",1u"mm"))
            b=MonoLightField2D([1.0 -2; 3im -5],wavelength=589u"nm",size=(1u"mm",1u"mm"))
            op(a,b)
        end
        @test_throws DimensionMismatch begin
            a=MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            b=MonoLightField2D([1.0 -2; 3im -5],wavelength=632.8u"nm",size=(1.1u"mm",1.1u"mm"))
            op(a,b)
        end
    end
    @testset "Multiplicative Operator" begin
        @test begin
            a=2+3im
            b=MonoLightField2D([1.0 -2; 3im -1+4im],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            c=MonoLightField2D([2.0+3.0im -4.0-6.0im; -9.0+6.0im -14.0+5.0im],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            a*b≈b*a≈c
        end
        @test begin
            a=MonoLightField2D([1.0 -2; 3im -1+4im],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            b=1.5+0.5im
            c=MonoLightField2D([0.6-0.2im -1.2+0.4im; 0.6+1.8im 0.2+2.6im],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            a/b≈b\a≈c
        end
    end
    @testset "Array-like" begin
        a=MonoLightField2D(rand(4,5),wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
        @test size(a) == size(a.data)
        @test a[1,1] == a.data[1,1]
        @test begin
            a[2,3] = 10
            a[2,3] == 10
        end
        @test abs.(a) == MonoLightField2D(a, data=abs.(a.data))
    end
    @testset "Plot Recipes" begin
        import RecipesBase.is_key_supported # We hack the RecipesBase to avoid MethodError when apply_recipe
        RecipesBase.is_key_supported(k::Symbol) = false
        l=MonoLightField2D(ones(4,5),wavelength=632.8u"nm",size=(1u"mm",1.2u"mm"))
        rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), l)
        @test getfield(rec[1],1)[:aspect_ratio] == :equal
        @test getfield(rec[1],1)[:size] == (4, 5)
        @test getfield(rec[1],1)[:xguide] == getfield(rec[1],1)[:yguide] == u"mm"
        @test rec[1].args == (-0.5:0.25:0.25, -0.6:0.24:0.36, ones(4,5))
    end
    @testset "Scale Invariant Propagation: $p" for p in [angularspectrum, fresnel2]
        @test begin # plane wave
            λ=632.8e-9
            lx=1.2e-3
            ly=1e-3
            d=1e-2
            Uin = ones(ComplexF64, 1024, 768)
            Uout = Uin*exp(2im*pi*d/λ)
            p(Uin, d, λ, lx, ly)≈Uout
        end
        @test begin # plane wave
            λ=632.8u"nm"
            lx=1.2u"mm"
            ly=1u"mm"
            d=1u"cm"
            Uin = MonoLightField2D(ones(ComplexF64, 1024, 768),wavelength=λ,size=(lx,ly))
            Uout = Uin*exp(2im*pi*convert(Float64, d/λ))
            p(Uin, d)≈Uout
        end
        @test begin # distance superposition
            λ=632.8u"nm"
            lx=ly=1u"mm"
            d1=1u"cm"
            d2=25u"mm"
            Uin = MonoLightField2D(rand(ComplexF64, 512, 512),wavelength=λ,size=(lx,ly))
            (p(d2)∘p(d1))(Uin)≈(p(d1)∘p(d2))(Uin)≈p(Uin, d1+d2)
        end
    end
    @testset "Fresnel S-FFT" begin
        @test begin # sample rate
            λ=632.8e-9
            lx=1.2e-3
            ly=1e-3
            d=1e-2
            nx = 16
            ny = 9
            Uin = ones(ComplexF64, nx, ny)
            (_, (lxo, lyo)) = fresnel1(Uin, d, λ, lx, ly)
            lxo≈λ*d*nx/lx && lyo≈λ*d*ny/ly
        end
    end
end
