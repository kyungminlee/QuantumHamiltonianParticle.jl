using QuantumHamiltonian
using LatticeTools
using QuantumHamiltonianParticle
using Test
using LinearAlgebra

@testset "Hubbard" begin
    electron_up = Fermion("↑")
    electron_down = Fermion("↓")

    particle_sector = ParticleSector(electron_up, electron_down)
    site = ParticleSite([
        ParticleState(particle_sector, "__", [0, 0], ( 0, 0)),
        ParticleState(particle_sector, "↑_", [1, 0], ( 1, 1)),
        ParticleState(particle_sector, "_↓", [0, 1], ( 1,-1)),
        ParticleState(particle_sector, "↑↓", [1, 1], ( 2, 0)),
    ])
    nsites = 4
    hilbert_space = ParticleHilbertSpace([site for i in 1:nsites])

    c_up_dag(i) = ParticleLadderUnit(particle_sector, 1, i, CREATION)
    c_up(i) = ParticleLadderUnit(particle_sector, 1, i, ANNIHILATION)

    c_dn_dag(i) = ParticleLadderUnit(particle_sector, 2, i, CREATION)
    c_dn(i) = ParticleLadderUnit(particle_sector, 2, i, ANNIHILATION)

    interaction_hamiltonian = sum(
        c_up_dag(i) * c_up(i) * c_dn_dag(i) * c_dn(i) for i in 1:nsites
    ) |> simplify

    hopping_hamiltonian = sum(
        let j = mod(i, nsites) + 1
            c_up_dag(i) * c_up(j) + c_dn_dag(i) * c_dn(j) + c_up_dag(j) * c_up(i) + c_dn_dag(j) * c_dn(i)
        end
            for i in 1:nsites
    ) |> simplify

    @testset "symmetry" begin
        t = SitePermutation([2,3,4,1])
        @test isinvariant(t, embed(hilbert_space, interaction_hamiltonian))
        @test isinvariant(t, embed(hilbert_space, hopping_hamiltonian))
        @test !isinvariant(t, embed(hilbert_space, c_up_dag(1)*c_up(1)))
    end

end
