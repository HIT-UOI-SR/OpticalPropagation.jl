using OpticalPropagation
using Test
using Unitful

@testset "OpticalPropagation.jl" begin
    @testset "MonoLightField2D Constructors" begin
        @test begin
            o=MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            ro=MonoLightField2D(o,data=[3 1; 2 4])
            r=MonoLightField2D([3 1; 2 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            ro≈r
        end
        @test begin
            o=MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            ro=MonoLightField2D(o,wavelength=1.15u"μm")
            r=MonoLightField2D([1 2; 3 4],wavelength=1.15u"μm",size=(1u"mm",1u"mm"))
            ro≈r
        end
        @test begin
            o=MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            ro=MonoLightField2D(o,size=(1.5u"mm",2.0u"mm"))
            r=MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1.5u"mm",2.0u"mm"))
            ro≈r
        end
        @test begin
            o=MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
            ro=MonoLightField2D(o,data=[3 1; 2 4],wavelength=1.15u"μm")
            r=MonoLightField2D([3 1; 2 4],wavelength=1.15u"μm",size=(1u"mm",1u"mm"))
            ro≈r
        end
        @test begin
            o=MonoLightField2D([1 2; 3 4],wavelength=632.8u"nm",size=(1u"mm",1u"mm"))
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
    @testset "Scale Invariant Propagation: $p" for p in [angularspectrum, fresnel2]
        @test begin
            λ=632.8e-9
            lx=ly=1e-3
            d=1e-2
            Uin = ones(ComplexF64, 1024, 1024)
            Uout = Uin*exp(2im*pi*d/λ)
            p(Uin, d, λ, lx, ly)≈Uout
        end
        @test begin
            λ=632.8u"nm"
            lx=ly=1u"mm"
            d=1u"cm"
            Uin = MonoLightField2D(ones(ComplexF64, 1024, 1024),wavelength=λ,size=(lx,ly))
            Uout = Uin*exp(2im*pi*convert(Float64, d/λ))
            p(Uin, d)≈Uout
        end
    end
end
