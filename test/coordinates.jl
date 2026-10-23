using FFTW

@testset "Physical coordinates" begin
    @testset "Legacy centred coordinates: $dims" for dims in ((4, 6), (5, 3), (1, 5), (3, 1), (1, 1))
        field = MonoLightField2D(ones(dims...); wavelength=633u"nm", size=(2u"mm", 3u"mm"))
        @test physicalsize(field) == field.size == (2u"mm", 3u"mm")
        @test gridcenter(field) == (0u"m", 0u"m")
        @test all(samplepitch(field) .≈ (2u"mm"/dims[1], 3u"mm"/dims[2]))
        @test axes(field) == axes(parent(field))
        @test parent(field) === field.data
        xs, ys = @inferred coordinates(field)
        @test xs isa AbstractRange
        @test ys isa AbstractRange
        @test length(xs) == dims[1]
        @test length(ys) == dims[2]
        @test collect(xs) ≈ [-1u"mm" + (i-1)*2u"mm"/dims[1] for i in 1:dims[1]]
        @test collect(ys) ≈ [-1.5u"mm" + (j-1)*3u"mm"/dims[2] for j in 1:dims[2]]
        @test first(xs) + length(xs)*step(xs) ≈ 1u"mm"
        @test first(ys) + length(ys)*step(ys) ≈ 1.5u"mm"
        @test @inferred(physicalsize(field)) == field.size
        @test @inferred(samplepitch(field)) == field.size ./ dims
        @test @inferred(gridcenter(field)) == field.center
        @test @inferred(coordinates(field, 1)) == xs
    end

    source = [complex(i, -2j) for i in 1:4, j in 1:3]
    field = MonoLightField2D(source; wavelength=633u"nm", size=(8u"mm", 6u"mm"),
                             center=(1u"cm", -3000u"μm"))
    @testset "Constructors and ownership" begin
        @test gridcenter(field) == (10u"mm", -3u"mm")
        @test collect(coordinates(field, 1)) ≈ [6, 8, 10, 12] .* u"mm"
        @test collect(coordinates(field, 2)) ≈ [-6, -4, -2] .* u"mm"
        @test parent(field) !== source
        snapshot = copy(parent(field))
        source[1, 1] = 1000
        @test parent(field) == snapshot
        for duplicate in (copy(field), MonoLightField2D(field))
            @test duplicate ≈ field
            @test parent(duplicate) !== parent(field)
            duplicate[1, 1] = 999
            @test parent(field) == snapshot
        end
        resized = MonoLightField2D(field; size=(4u"mm", 3u"mm"))
        @test gridcenter(resized) == gridcenter(field)
        @test samplepitch(resized) == samplepitch(field) ./ 2
        resampled = MonoLightField2D(field; data=ones(8, 6))
        @test physicalsize(resampled) == physicalsize(field)
        @test gridcenter(resampled) == gridcenter(field)
        @test samplepitch(resampled) == samplepitch(field) ./ 2
        moved = MonoLightField2D(field; center=(0u"m", 0u"m"))
        @test gridcenter(moved) == (0u"m", 0u"m")
        @test !(moved ≈ field)
        @test occursin("center=", repr(field))
        @test !occursin("center=", repr(moved))
        @test occursin("window centre:", sprint(show, MIME"text/plain"(), field))
        @test_throws ArgumentError coordinates(field, 0)
        @test_throws ArgumentError coordinates(field, 3)
        @test_throws ArgumentError MonoLightField2D(field; center=(NaN*u"m", 0u"m"))
        @test_throws ArgumentError MonoLightField2D(field; center=(0u"m", Inf*u"m"))
        emptyfield = MonoLightField2D(zeros(0, 3); wavelength=633u"nm", size=(1u"mm", 1u"mm"))
        @test_throws ArgumentError samplepitch(emptyfield)
        @test_throws ArgumentError coordinates(emptyfield, 1)
    end

    @testset "Arithmetic, fused broadcasting and grid checks" begin
        for result in (2field, field*2, field/2, 2\field, abs.(field), field .+ 1,
                       2 .* field .+ field .* ones(size(field)))
            @test gridcenter(result) == gridcenter(field)
            @test physicalsize(result) == physicalsize(field)
            @test result.wavelength == field.wavelength
        end
        @test parent(big"1.5"*field) == ComplexF64.(big"1.5" .* parent(field))
        @test parent(field + field) == 2 .* parent(field)
        @test parent(field - field) == zeros(size(field))
        @test parent(2 .* field .+ field .* ones(size(field))) == 3 .* parent(field)
        @test parent(field .+ [1, 2, 3, 4]) == parent(field) .+ [1, 2, 3, 4]
        @test isapprox(copy(field), field; rtol=1e-12)
        # A labelled grid must not acquire extra samples through broadcasting.
        thin = MonoLightField2D(field; data=ones(1, 3))
        @test_throws DimensionMismatch thin .+ ones(4, 3)
        @test_throws DimensionMismatch field .+ ones(4, 3, 2)

        incompatible = (
            MonoLightField2D(field; center=(11u"mm", -3u"mm")),
            MonoLightField2D(field; size=(9u"mm", 6u"mm")),
            MonoLightField2D(field; wavelength=532u"nm"),
            thin,
        )
        for other in incompatible
            @test_throws DimensionMismatch field + other
            @test_throws DimensionMismatch field - other
            @test_throws DimensionMismatch field .+ other
            @test_throws DimensionMismatch 2 .* field .+ abs.(other)
            dest = copy(field)
            original = copy(parent(dest))
            @test_throws DimensionMismatch dest .= other
            @test parent(dest) == original
            @test_throws DimensionMismatch dest .= 2 .* field .+ other
            @test parent(dest) == original
        end
        dest = copy(field)
        dest .= 2 .* dest .+ field
        @test parent(dest) == 3 .* parent(field)
        @test gridcenter(dest) == gridcenter(field)
        dest .= 4
        @test all(==(4), parent(dest))
        dest .= ones(size(dest)) .+ 2
        @test all(==(3), parent(dest))
        dest .= field
        @test dest ≈ field
        # Check all labelled inputs even when the destination is a plain matrix.
        raw = zeros(ComplexF64, size(field))
        raw .= field .+ field
        @test raw == 2 .* parent(field)
        @test_throws DimensionMismatch raw .= field .+ first(incompatible)

        far = MonoLightField2D(field; center=(1e9u"m", 0u"m"))
        shifted = MonoLightField2D(far; center=(1e9u"m" + samplepitch(far)[1], 0u"m"))
        @test_throws DimensionMismatch far + shifted
        @test !(far ≈ shifted)
        # Equivalent units and roundoff-sized displacements remain compatible.
        equivalent = MonoLightField2D(field; center=(10000u"μm", -0.3u"cm"))
        @test field + equivalent ≈ 2field
    end

    @testset "Coordinate-preserving crops" begin
        result = crop(field, 2:3, 2:3)
        @test parent(result) == parent(field)[2:3, 2:3]
        @test samplepitch(result) == samplepitch(field)
        @test physicalsize(result) == (4u"mm", 4u"mm")
        @test gridcenter(result) == (10u"mm", -2u"mm")
        @test coordinates(result, 1) ≈ coordinates(field, 1)[2:3]
        @test coordinates(result, 2) ≈ coordinates(field, 2)[2:3]
        @test result.wavelength == field.wavelength
        original = copy(parent(field))
        result[1, 1] = 999
        @test parent(field) == original
        @test crop(field, :, :) ≈ field
        @test parent(crop(field, :, :)) !== parent(field)
        @test coordinates(crop(field, :, 1:1), 2) ≈ coordinates(field, 2)[1:1]
        @test field[2:3, 2:3] isa Matrix
        @test_throws BoundsError crop(field, 0:2, 1:2)
        @test_throws BoundsError crop(field, 1:5, 1:2)
        @test_throws ArgumentError crop(field, 2:1, :)
        @test_throws MethodError crop(field, 1:2:3, :)
    end

    @testset "Plot coordinates" begin
        rec = RecipesBase.apply_recipe(Dict{Symbol,Any}(), field)[1]
        @test rec.args[1] ≈ [6, 8, 10, 12]
        @test rec.args[2] ≈ [-6, -4, -2]
        @test rec.args[3] == abs.(parent(field))
        @test rec.plotattributes[:aspect_ratio] == :equal
        metre_recipe = RecipesBase.apply_recipe(Dict{Symbol,Any}(:unit => u"m"), field)[1]
        @test metre_recipe.args[1] ≈ [0.006, 0.008, 0.010, 0.012]
        @test metre_recipe.args[2] ≈ [-0.006, -0.004, -0.002]
    end

    @testset "Matrix FFT access and coordinate storage" begin
        @test fft(parent(field)) ≈ fft(copy(parent(field)))
        @test parent(field) === field.data
        small = MonoLightField2D(ones(2, 3); wavelength=633u"nm", size=(1u"mm", 1u"mm"))
        large = MonoLightField2D(ones(128, 129); wavelength=633u"nm", size=small.size)
        @test Base.summarysize(coordinates(small)) == Base.summarysize(coordinates(large))
    end

    @testset "Propagation metadata without numerical kernel changes" begin
        for dims in ((4, 6), (5, 3))
            input = MonoLightField2D([complex(i, -j) for i in 1:dims[1], j in 1:dims[2]];
                wavelength=633u"nm", size=(0.1u"mm", 0.08u"mm"), center=(1u"mm", -2u"mm"))
            centred = MonoLightField2D(input; center=(0u"m", 0u"m"))
            for propagate in (angularspectrum, fresnel2, fresnelkirchhoff)
                output = propagate(input, 2u"mm")
                @test parent(output) == parent(propagate(centred, 2u"mm"))
                @test gridcenter(output) == gridcenter(input)
                @test physicalsize(output) == physicalsize(input)
                @test size(output) == size(input)
                @test output.wavelength == input.wavelength
                @test propagate(2u"mm")(input) ≈ output
            end
            for distance in (2u"mm", -2u"mm")
                output = fresnel1(input, distance)
                raw, widths = fresnel1(parent(input), OpticalPropagation.auval(distance),
                    OpticalPropagation.auval(input.wavelength),
                    OpticalPropagation.auval(input.size[1]), OpticalPropagation.auval(input.size[2]))
                @test parent(output) == raw
                @test physicalsize(output) == widths .* u"m"
                @test gridcenter(output) == gridcenter(input)
                @test output.wavelength == input.wavelength
                @test fresnel1(distance)(input) ≈ output
            end
            identity = fresnelkirchhoff(input, 0u"mm")
            @test identity ≈ input
            @test parent(identity) !== parent(input)
            @test_throws ArgumentError fresnelkirchhoff(input, 2u"mm"; radius=10u"cm")
            @test_throws ArgumentError fresnelkirchhoff(input, 2u"mm"; radius=-10u"cm")
            @test_throws ArgumentError fresnelkirchhoff(2u"mm"; radius=10u"cm")(input)
        end
    end
end
