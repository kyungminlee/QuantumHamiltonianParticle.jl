using Test
using LatticeTools
using ExactDiagonalization
using Particle

@testset "hilbert" begin
    p = ParticleSector(Boson(:m, 2), Fermion(:f))

    @testset "exception" begin
        states = [
            ParticleState(p, "↑.", [0, 0], ( 1, 0), UInt8),
            ParticleState(p, "0.", [1, 0], ( 0, 0), UInt8),
            ParticleState(p, "↑.", [2, 0], (-1, 0), UInt8),
            ParticleState(p, "↑f", [0, 1], ( 1, 1), UInt8),
            ParticleState(p, "0f", [1, 1], ( 0, 1), UInt8),
            ParticleState(p, "↑f", [2, 1], (-1, 1), UInt8),
        ]
        site = ParticleSite(states)
        ParticleHilbertSpace([site, site])
        @test_throws ArgumentError ParticleHilbertSpace([site, site, site])
    end

    states = [
        ParticleState(p, "↑.", [0, 0], ( 1, 0)),
        ParticleState(p, "0.", [1, 0], ( 0, 0)),
        ParticleState(p, "↑.", [2, 0], (-1, 0)),
        ParticleState(p, "↑f", [0, 1], ( 1, 1)),
        ParticleState(p, "0f", [1, 1], ( 0, 1)),
        ParticleState(p, "↑f", [2, 1], (-1, 1)),
    ]
    site = ParticleSite(states)
    hilbert = ParticleHilbertSpace([site, site, site, site])
    @test hilbert == ParticleHilbertSpace([site, site, site, site])
    @test hilbert != ParticleHilbertSpace([site, site, site])

    @test scalartype(hilbert) == Bool
    @test valtype(hilbert) == Bool
    @test qntype(hilbert) == Tuple{Int, Int}
    @test basespace(hilbert) == hilbert

    @testset "particlesector" begin
        @test exchangesign(hilbert, 1) == 1
        @test exchangesign(hilbert, 2) == -1
        @test numspecies(hilbert) == 2
        @test speciescount(hilbert) == 2
        @test getspecies(hilbert, 1) == Boson{:m, 2}
        @test getspecies(hilbert, 2) == Fermion{:f}
        @test getspeciesname(hilbert, 1) == :m
        @test getspeciesname(hilbert, 2) == :f
    end


    @test bitwidth(hilbert) == 12

    @test bitoffset(hilbert, 1, 1) == 0
    @test bitoffset(hilbert, 2, 1) == 2
    @test bitoffset(hilbert, 1, 2) == 3
    @test bitoffset(hilbert, 2, 2) == 5
    @test bitoffset(hilbert, 1) == 0
    @test bitoffset(hilbert, 2) == 3

    @test get_bitmask(hilbert, 1, 1) == 0b000000000011
    @test get_bitmask(hilbert, 2, 1) == 0b000000000100
    @test get_bitmask(hilbert, 1, 2) == 0b000000011000
    @test get_bitmask(hilbert, 2, 2) == 0b000000100000

    @test get_bitmask(hilbert, 1, :) == 0b011011011011
    @test get_bitmask(hilbert, :, :) == 0b111111111111
    @test get_bitmask(hilbert, :, 2) == 0b000000111000
    # @test get_bitmask(hilbert, :, [2,3]) == 0b000111111000
    # TODO(kyungminlee 20200917) tests for colon and vector

    @test get_bitmask(hilbert, 1, [2,3])    == 0b000011011000
    @test get_bitmask(hilbert, [1], 2)      == 0b000000011000
    @test get_bitmask(hilbert, [1,2], 2)    == 0b000000111000
    @test get_bitmask(hilbert, [1,], [2,3]) == 0b000011011000

    @test_throws BoundsError get_bitmask(hilbert, 3, 1)
    @test_throws BoundsError get_bitmask(hilbert, 1, 8)

    @test_throws BoundsError get_bitmask(hilbert, 3, :)
    @test_throws BoundsError get_bitmask(hilbert, :, 8)

    @test_throws BoundsError get_bitmask(hilbert, 3, [1,2])
    @test_throws BoundsError get_bitmask(hilbert, [1,2], 8)




    # @test get_bitmask(hi)
    # @show hilbert
end
