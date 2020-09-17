using Test
using Particle
using LatticeTools
using ExactDiagonalization


@testset "state-site" begin
    p = ParticleSector(Boson(:m, 2), Fermion(:f))
    @testset "state" begin
        state = ParticleState(p, "↑f", [0, 1], ( 1, 1))
        @test state.name == "↑f"
        @test state.occupancy_binary == UInt(0b100)
        @test state.quantum_number == (1, 1)

        state1 = ParticleState(p, "↑f", [0, 1], ( 1, 1))
        state2 = ParticleState(p, "↓f", [1, 1], (-1, 1))

        @test state == state1
        @test state != state2

        @test qntype(state) == Tuple{Int, Int}

        @test speciescount(state) == 2
        @test numspecies(state) == 2

    end
end
