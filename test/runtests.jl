using OpticalPropagation
using Test
using Unitful

@testset "OpticalPropagation.jl" begin
    @testset "MonoLightField2D" begin
        @test begin
            o=MonoLightField2D([1 2; 3 4],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            ro=MonoLightField2D(o,distribution_data=[3 1; 2 4])
            r=MonoLightField2D([3 1; 2 4],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            ro≈r
        end
        @test begin
            o=MonoLightField2D([1 2; 3 4],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            ro=MonoLightField2D(o,wavelength=1.15*u"μm")
            r=MonoLightField2D([1 2; 3 4],wavelength=1.15*u"μm",size=(1*u"mm",1*u"mm"))
            ro≈r
        end
        @test begin
            o=MonoLightField2D([1 2; 3 4],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            ro=MonoLightField2D(o,size=(1.5*u"mm",2.0*u"mm"))
            r=MonoLightField2D([1 2; 3 4],wavelength=632.8*u"nm",size=(1.5*u"mm",2.0*u"mm"))
            ro≈r
        end
        @test begin
            o=MonoLightField2D([1 2; 3 4],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            ro=MonoLightField2D(o,distribution_data=[3 1; 2 4],wavelength=1.15*u"μm")
            r=MonoLightField2D([3 1; 2 4],wavelength=1.15*u"μm",size=(1*u"mm",1*u"mm"))
            ro≈r
        end
        @test begin
            o=MonoLightField2D([1 2; 3 4],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            ro=MonoLightField2D(o,wavelength=1.15*u"μm",size=(1.5*u"mm",2.0*u"mm"))
            r=MonoLightField2D([1 2; 3 4],wavelength=1.15*u"μm",size=(1.5*u"mm",2.0*u"mm"))
            ro≈r
        end
        @test begin
            a=MonoLightField2D([1 2; 3 4],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            b=MonoLightField2D([1.0 -2; 3im -5],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            c=MonoLightField2D([2.0 0.0; 3.0+3.0im -1.0],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            a+b≈c
        end
        @test begin
            a=MonoLightField2D([1 2; 3 4],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            b=MonoLightField2D([1.0 -2; 3im -5],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            c=MonoLightField2D([0.0 4.0; 3.0-3.0im 9.0],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            a-b≈c
        end
        @test begin
            a=2+3im
            b=MonoLightField2D([1.0 -2; 3im -1+4im],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            c=MonoLightField2D([2.0+3.0im -4.0-6.0im; -9.0+6.0im -14.0+5.0im],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            a*b≈b*a≈c
        end
        @test begin
            a=MonoLightField2D([1.0 -2; 3im -1+4im],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            b=1.5+0.5im
            c=MonoLightField2D([0.6-0.2im -1.2+0.4im; 0.6+1.8im 0.2+2.6im],wavelength=632.8*u"nm",size=(1*u"mm",1*u"mm"))
            a/b≈b\a≈c
        end
    end
end
