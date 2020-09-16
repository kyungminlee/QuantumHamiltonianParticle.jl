using ExactDiagonalization
using Particle
using Test
using LinearAlgebra
# using Formatting

@testset "Hubbard" begin
    electron_up = Fermion{Symbol("↑")}()
    electron_down = Fermion{Symbol("↓")}()

    particle_sector = make_particle_sector(electron_up, electron_down)
    site = ParticleSite([
        ParticleState(particle_sector, "__", [0, 0], ( 0, 0)),
        ParticleState(particle_sector, "↑_", [1, 0], ( 1, 1)),
        ParticleState(particle_sector, "_↓", [0, 1], ( 1,-1)),
        ParticleState(particle_sector, "↑↓", [1, 1], ( 2, 0)),
    ])
    nsites = 4
    hs = ParticleHilbertSpace([site for i in 1:nsites])

    c_up_dag(i) = LadderUnitOperator(particle_sector, 1, i, CREATION)
    c_up(i) = LadderUnitOperator(particle_sector, 1, i, ANNIHILATION)

    c_dn_dag(i) = LadderUnitOperator(particle_sector, 2, i, CREATION)
    c_dn(i) = LadderUnitOperator(particle_sector, 2, i, ANNIHILATION)


    interaction_hamiltonian = sum(c_up_dag(i) * c_up(i) * c_dn_dag(i) * c_dn(i) for i in 1:nsites) |> simplify
    hopping_hamiltonian = sum(
        let j = mod(i, nsites) + 1
            c_up_dag(i) * c_up(j) + c_dn_dag(i) * c_dn(j) + c_up_dag(j) * c_up(i) + c_dn_dag(j) * c_dn(i)
        end
            for i in 1:nsites
    ) |> simplify

    # prettyprintln( interaction_hamiltonian )
    # prettyprintln( hopping_hamiltonian )

end
