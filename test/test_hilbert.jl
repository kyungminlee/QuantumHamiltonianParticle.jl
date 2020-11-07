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
    hilbert = ParticleHilbertSpace([site, site, site])
    @test hilbert != ParticleHilbertSpace([site, site, site, site])
    @test hilbert == ParticleHilbertSpace([site, site, site])

    for h in [hilbert, typeof(hilbert)]
        @test scalartype(h) == Bool
        @test valtype(h) == Bool
        @test qntype(h) == Tuple{Int, Int}
    end

    @test basespace(hilbert) == hilbert

    @testset "particlesector" begin
        for h in [hilbert, typeof(hilbert)]
            @test exchangesign(h, 1) == 1
            @test exchangesign(h, 2) == -1
            @test numspecies(h) == 2
            @test speciescount(h) == 2
            @test getspecies(h, 1) == Boson{:m, 2}
            @test getspecies(h, 2) == Fermion{:f}
            @test getspeciesname(h, 1) == :m
            @test getspeciesname(h, 2) == :f
        end
    end

    @test bitwidth(hilbert) == 9

    @test bitoffset(hilbert, 1, 1) == 0
    @test bitoffset(hilbert, 2, 1) == 2
    @test bitoffset(hilbert, 1, 2) == 3
    @test bitoffset(hilbert, 2, 2) == 5
    @test bitoffset(hilbert, 1) == 0
    @test bitoffset(hilbert, 2) == 3

    @test get_bitmask(hilbert, 1, 1) == 0b000000011
    @test get_bitmask(hilbert, 2, 1) == 0b000000100
    @test get_bitmask(hilbert, 1, 2) == 0b000011000
    @test get_bitmask(hilbert, 2, 2) == 0b000100000

    @test get_bitmask(hilbert, 1, :) == 0b011011011
    @test get_bitmask(hilbert, :, :) == 0b111111111
    @test get_bitmask(hilbert, :, 2) == 0b000111000
    # @test get_bitmask(hilbert, :, [2,3]) == 0b000111111000
    # TODO(kyungminlee 20200917) tests for colon and vector

    @test get_bitmask(hilbert, 1, [2,3])    == 0b011011000
    @test get_bitmask(hilbert, [1], 2)      == 0b000011000
    @test get_bitmask(hilbert, [1,2], 2)    == 0b000111000
    @test get_bitmask(hilbert, [1,], [2,3]) == 0b011011000

    @test get_bitmask(hilbert) == 0b111111111

    @test_throws BoundsError get_bitmask(hilbert, 3, 1)
    @test_throws BoundsError get_bitmask(hilbert, 1, 8)

    @test_throws BoundsError get_bitmask(hilbert, 3, :)
    @test_throws BoundsError get_bitmask(hilbert, :, 8)

    @test_throws BoundsError get_bitmask(hilbert, 3, [1,2])
    @test_throws BoundsError get_bitmask(hilbert, [1,2], 8)

    @test get_parity_bitmask(hilbert, 2, 3) == 0b000100100  # parity bits for fermions
    @test get_parity_bitmask(hilbert, 1, 3) == 0b000000000  # nothing for bosons

    @test get_occupancy(hilbert, 1, 1, 0b001_110_100) == 0
    @test get_occupancy(hilbert, 2, 1, 0b001_110_100) == 1
    @test get_occupancy(hilbert, 1, 2, 0b001_110_100) == 2
    @test get_occupancy(hilbert, 2, 2, 0b001_110_100) == 1
    @test get_occupancy(hilbert, 1, 3, 0b001_110_100) == 1
    @test get_occupancy(hilbert, 2, 3, 0b001_110_100) == 0

    @test set_occupancy(hilbert, 1, 3, 0b001_110_100, 0) == 0b000_110_100
    @test_throws ArgumentError set_occupancy(hilbert, 1, 3, 0b001_110_100, 3)

    # three sites
    @test compress(hilbert, CartesianIndex(2, 1, 6)) == UInt(0b110_000_001)
    @test extract(hilbert, 0b110_000_001) == CartesianIndex(2, 1, 6)
    @test hilbert[CartesianIndex(2, 1, 6)] == UInt(0b110_000_001)
    @test hilbert[2, 1, 6] == UInt(0b110_000_001)

    #=
    @test get_bitmask(hilbert, 1, 1) == 0b000000011
    @test get_bitmask(hilbert, 2, 1) == 0b000000100
    @test get_bitmask(hilbert, 1, 2) == 0b000011000
    @test get_bitmask(hilbert, 2, 2) == 0b000100000
    =#

    @test Set(quantum_number_sectors(hilbert)) == Set([(sz, c) for sz in -3:3, c in 0:3])
    @test get_quantum_number(hilbert, 0b001_110_100) == (0-1+1, 1+1)
    @test get_quantum_number(hilbert, [4, 6, 2]) == (0-1+1, 1+1)

    @test keys(hilbert) == CartesianIndices((1:6, 1:6, 1:6))

    basis_list = ExactDiagonalization.hs_get_basis_list(hilbert)
    @test basis_list == [
        0b000000000, 0b000000001, 0b000000010, 0b000000100, 0b000000101, 0b000000110,
        0b000001000, 0b000001001, 0b000001010, 0b000001100, 0b000001101, 0b000001110,
        0b000010000, 0b000010001, 0b000010010, 0b000010100, 0b000010101, 0b000010110,
        0b000100000, 0b000100001, 0b000100010, 0b000100100, 0b000100101, 0b000100110,
        0b000101000, 0b000101001, 0b000101010, 0b000101100, 0b000101101, 0b000101110,
        0b000110000, 0b000110001, 0b000110010, 0b000110100, 0b000110101, 0b000110110,
        0b001000000, 0b001000001, 0b001000010, 0b001000100, 0b001000101, 0b001000110,
        0b001001000, 0b001001001, 0b001001010, 0b001001100, 0b001001101, 0b001001110,
        0b001010000, 0b001010001, 0b001010010, 0b001010100, 0b001010101, 0b001010110,
        0b001100000, 0b001100001, 0b001100010, 0b001100100, 0b001100101, 0b001100110,
        0b001101000, 0b001101001, 0b001101010, 0b001101100, 0b001101101, 0b001101110,
        0b001110000, 0b001110001, 0b001110010, 0b001110100, 0b001110101, 0b001110110,
        0b010000000, 0b010000001, 0b010000010, 0b010000100, 0b010000101, 0b010000110,
        0b010001000, 0b010001001, 0b010001010, 0b010001100, 0b010001101, 0b010001110,
        0b010010000, 0b010010001, 0b010010010, 0b010010100, 0b010010101, 0b010010110,
        0b010100000, 0b010100001, 0b010100010, 0b010100100, 0b010100101, 0b010100110,
        0b010101000, 0b010101001, 0b010101010, 0b010101100, 0b010101101, 0b010101110,
        0b010110000, 0b010110001, 0b010110010, 0b010110100, 0b010110101, 0b010110110,
        0b100000000, 0b100000001, 0b100000010, 0b100000100, 0b100000101, 0b100000110,
        0b100001000, 0b100001001, 0b100001010, 0b100001100, 0b100001101, 0b100001110,
        0b100010000, 0b100010001, 0b100010010, 0b100010100, 0b100010101, 0b100010110,
        0b100100000, 0b100100001, 0b100100010, 0b100100100, 0b100100101, 0b100100110,
        0b100101000, 0b100101001, 0b100101010, 0b100101100, 0b100101101, 0b100101110,
        0b100110000, 0b100110001, 0b100110010, 0b100110100, 0b100110101, 0b100110110,
        0b101000000, 0b101000001, 0b101000010, 0b101000100, 0b101000101, 0b101000110,
        0b101001000, 0b101001001, 0b101001010, 0b101001100, 0b101001101, 0b101001110,
        0b101010000, 0b101010001, 0b101010010, 0b101010100, 0b101010101, 0b101010110,
        0b101100000, 0b101100001, 0b101100010, 0b101100100, 0b101100101, 0b101100110,
        0b101101000, 0b101101001, 0b101101010, 0b101101100, 0b101101101, 0b101101110,
        0b101110000, 0b101110001, 0b101110010, 0b101110100, 0b101110101, 0b101110110,
        0b110000000, 0b110000001, 0b110000010, 0b110000100, 0b110000101, 0b110000110,
        0b110001000, 0b110001001, 0b110001010, 0b110001100, 0b110001101, 0b110001110,
        0b110010000, 0b110010001, 0b110010010, 0b110010100, 0b110010101, 0b110010110,
        0b110100000, 0b110100001, 0b110100010, 0b110100100, 0b110100101, 0b110100110,
        0b110101000, 0b110101001, 0b110101010, 0b110101100, 0b110101101, 0b110101110,
        0b110110000, 0b110110001, 0b110110010, 0b110110100, 0b110110101, 0b110110110,
    ]
    @test length(basis_list) == 6^3
    @test_throws ArgumentError ExactDiagonalization.hs_get_basis_list(hilbert, UInt8)
end
