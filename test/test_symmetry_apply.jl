using Test
using Particle
using LatticeTools
using ExactDiagonalization

@testset "symmetry" begin
    perm = SitePermutation([2,3,1])  # 1->2, 2->3, 3->1
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

    @test symmetry_apply(perm, c(1,1)) == c(1,2)
    @test symmetry_apply(perm, c(1,3)) == c(1,1)
    @test symmetry_apply(perm, c(2,3)) == c(2,1)
    @test symmetry_apply(perm, cdag(1,2)) == cdag(1,3)
    @test symmetry_apply(perm, cdag(2,3)) == cdag(2,1)

    @test symmetry_apply(particle_hilbert_space, perm, c(1,1)) == c(1,2)
    @test symmetry_apply(particle_hilbert_space, perm, c(1,3)) == c(1,1)
    @test symmetry_apply(particle_hilbert_space, perm, c(2,3)) == c(2,1)
    @test symmetry_apply(particle_hilbert_space, perm, cdag(1,2)) == cdag(1,3)
    @test symmetry_apply(particle_hilbert_space, perm, cdag(2,3)) == cdag(2,1)

    @test symmetry_apply(perm, ec(1,1)) == ec(1,2)
    @test symmetry_apply(perm, ec(1,3)) == ec(1,1)
    @test symmetry_apply(perm, ec(2,3)) == ec(2,1)
    @test symmetry_apply(perm, ecdag(1,2)) == ecdag(1,3)
    @test symmetry_apply(perm, ecdag(2,3)) == ecdag(2,1)


    @test symmetry_apply(particle_hilbert_space, perm, cdag(1,2)*c(2,3)) == cdag(1,3)*c(2,1)
    @test symmetry_apply(particle_hilbert_space, perm, cdag(1,2)*c(2,3)+c(2,1)) == cdag(1,3)*c(2,1)+c(2,2)

    @test symmetry_apply(perm, cdag(1,2)*c(2,3)) == cdag(1,3)*c(2,1)
    @test symmetry_apply(perm, cdag(1,2)*c(2,3)+c(2,1)) == cdag(1,3)*c(2,1)+c(2,2)

    @test !isinvariant(perm, c(1,1) + c(1,2))
    @test isinvariant(perm, c(1,1) + c(1,2) + c(1,3))
end
