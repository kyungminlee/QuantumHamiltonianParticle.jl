using Test
using Particle
using LatticeTools
using ExactDiagonalization

@testset "symmetry" begin
    permutation = SitePermutation([2,3,1])  # 1->2, 2->3, 3->1
    particle_sector = ParticleSector(Fermion("f"), Boson("S=1", 2))
    particle_site = ParticleSite([
        ParticleState(particle_sector, "_↑", [0,0]),
        ParticleState(particle_sector, "f↑", [1,0]),
        ParticleState(particle_sector, "_⋅", [0,1]),
        ParticleState(particle_sector, "f⋅", [1,1]),
        ParticleState(particle_sector, "_↓", [0,2]),
        ParticleState(particle_sector, "f↓", [1,2]),
    ])
    particle_hilbert_space = ParticleHilbertSpace([particle_site, particle_site, particle_site])

    c(i,j) = ParticleLadderUnit(particle_sector, i, j, ANNIHILATION)
    cdag(i,j) = ParticleLadderUnit(particle_sector, i, j, CREATION)

    ec(i,j) = embed(particle_hilbert_space, ParticleLadderUnit(particle_sector, i, j, ANNIHILATION))
    ecdag(i,j) = embed(particle_hilbert_space, ParticleLadderUnit(particle_sector, i, j, CREATION))

    @testset "ladder operator" begin
        @test symmetry_apply(permutation, c(1,1)) == c(1,2)
        @test symmetry_apply(permutation, c(1,3)) == c(1,1)
        @test symmetry_apply(permutation, c(2,3)) == c(2,1)
        @test symmetry_apply(permutation, cdag(1,2)) == cdag(1,3)
        @test symmetry_apply(permutation, cdag(2,3)) == cdag(2,1)

        @test symmetry_apply(particle_hilbert_space, permutation, c(1,1)) == c(1,2)
        @test symmetry_apply(particle_hilbert_space, permutation, c(1,3)) == c(1,1)
        @test symmetry_apply(particle_hilbert_space, permutation, c(2,3)) == c(2,1)
        @test symmetry_apply(particle_hilbert_space, permutation, cdag(1,2)) == cdag(1,3)
        @test symmetry_apply(particle_hilbert_space, permutation, cdag(2,3)) == cdag(2,1)

        @test symmetry_apply(permutation, ec(1,1)) == ec(1,2)
        @test symmetry_apply(permutation, ec(1,3)) == ec(1,1)
        @test symmetry_apply(permutation, ec(2,3)) == ec(2,1)
        @test symmetry_apply(permutation, ecdag(1,2)) == ecdag(1,3)
        @test symmetry_apply(permutation, ecdag(2,3)) == ecdag(2,1)


        @test symmetry_apply(particle_hilbert_space, permutation, cdag(1,2)*c(2,3)) == cdag(1,3)*c(2,1)
        @test symmetry_apply(particle_hilbert_space, permutation, cdag(1,2)*c(2,3)+c(2,1)) == cdag(1,3)*c(2,1)+c(2,2)

        @test symmetry_apply(permutation, cdag(1,2)*c(2,3)) == cdag(1,3)*c(2,1)
        @test symmetry_apply(permutation, cdag(1,2)*c(2,3)+c(2,1)) == cdag(1,3)*c(2,1)+c(2,2)

        @test !isinvariant(permutation, c(1,1) + c(1,2))
        @test isinvariant(permutation, c(1,1) + c(1,2) + c(1,3))
    end

    @testset "state" begin
        @test symmetry_apply(particle_hilbert_space, permutation, 0b000_000_000) == (0b000_000_000, 1)

        @test symmetry_apply(particle_hilbert_space, permutation, 0b000_000_001) == (0b000_001_000, 1)
        @test symmetry_apply(particle_hilbert_space, permutation, 0b000_001_000) == (0b001_000_000, 1)
        @test symmetry_apply(particle_hilbert_space, permutation, 0b001_000_000) == (0b000_000_001, 1)

        @test symmetry_apply(particle_hilbert_space, permutation, 0b000_000_011) == (0b000_011_000, 1)
        @test symmetry_apply(particle_hilbert_space, permutation, 0b000_011_000) == (0b011_000_000, 1)
        @test symmetry_apply(particle_hilbert_space, permutation, 0b011_000_000) == (0b000_000_011, 1)

        @test symmetry_apply(particle_hilbert_space, permutation, 0b010_001_011) == (0b001_011_010, 1)
        @test symmetry_apply(particle_hilbert_space, permutation, 0b001_011_010) == (0b011_010_001, -1)
        @test symmetry_apply(particle_hilbert_space, permutation, 0b011_010_001) == (0b010_001_011, -1)
    end
end
