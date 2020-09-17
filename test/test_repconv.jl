using Test
using LatticeTools
using Particle
using ExactDiagonalization

@testset "repconv" begin
    p = ParticleSector(Boson(:m, 2), Fermion(:f))

    site = ParticleSite([
        ParticleState(p, "↑.", [0, 0], ( 1, 0)),
        ParticleState(p, "0.", [1, 0], ( 0, 0)),
        ParticleState(p, "↑.", [2, 0], (-1, 0)),
        ParticleState(p, "↑f", [0, 1], ( 1, 1)),
        ParticleState(p, "0f", [1, 1], ( 0, 1)),
        ParticleState(p, "↑f", [2, 1], (-1, 1)),
    ])

    @testset "site" begin
        @test state2occbin(site, 3) == 0b10
        #   state       occvec      occbin
        #     - occvec    - state     - state
        #     - occbin    - occbin    - occvec
        for (istate, state) in enumerate(site.states)
            ov = state2occvec(site, istate)
            ob = state2occbin(site, istate)

            istate1 = occvec2state(site, ov)
            ob1 = occvec2occbin(site, ov)
            @test istate1 == istate
            @test ob1 == ob

            istate2 = occbin2state(site, ob)
            ov2 = occbin2occvec(site, ob)
            @test istate2 == istate
            @test ov2 == ov
        end
    end

    @testset "hilbert" begin

    end
end
