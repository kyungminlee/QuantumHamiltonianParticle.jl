using Test
using LinearAlgebra
using ExactDiagonalization
using Particle

@testset "Ladder Iterator" begin
    p = ParticleSector(Boson(:m, 2), Fermion(:f))
    cdag(i, j) = ParticleLadderUnit(p, i, j, CREATION)
    c(i, j)    = ParticleLadderUnit(p, i, j, ANNIHILATION)
    site = ParticleSite([
        ParticleState(p, "__", [0, 0], (0, 0)),
        ParticleState(p, "b_", [1, 0], (1, 0)),
        ParticleState(p, "B_", [2, 0], (2, 0)),
        ParticleState(p, "_f", [0, 1], (0, 1)),
        ParticleState(p, "bf", [1, 1], (1, 1)),
        ParticleState(p, "Bf", [2, 1], (2, 1)),
    ])
    hs = ParticleHilbertSpace([site, site, site])
    @test collect(get_column_iterator(hs, cdag(1,2), UInt(0b000_000_000))) == [UInt(0b000_001_000) => 1.0]
    @test collect(get_column_iterator(hs, cdag(1,2), UInt(0b000_000_001))) == [UInt(0b000_001_001) => 1.0]
    @test collect(get_column_iterator(hs, cdag(2,2), UInt(0b000_000_000))) == [UInt(0b000_100_000) => 1.0]
    @test collect(get_column_iterator(hs, cdag(2,2), UInt(0b000_000_001))) == [UInt(0b000_100_001) => 1.0]
    @test collect(get_column_iterator(hs, cdag(2,2), UInt(0b000_000_101))) == [UInt(0b000_100_101) => -1.0]
end
