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


    @test bitwidth(hilbert) == 12

    @test bitoffset(hilbert, 1, 1) == 0
    @test bitoffset(hilbert, 2, 1) == 2
    @test bitoffset(hilbert, 1, 2) == 3
    @test bitoffset(hilbert, 2, 2) == 5
    @test bitoffset(hilbert, 1) == 0
    @test bitoffset(hilbert, 2) == 3

    @test get_bitmask(hilbert, 1, 1) == 0b000011
    @test get_bitmask(hilbert, 2, 1) == 0b000100
    @test get_bitmask(hilbert, 1, 2) == 0b011000
    @test get_bitmask(hilbert, 2, 2) == 0b100000

    # @test get_bitmask(hi)
    # @show hilbert
end
