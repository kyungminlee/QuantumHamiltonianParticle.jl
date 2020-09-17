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
        hilbert = ParticleHilbertSpace([site, site])
        for sv in keys(hilbert)
            sv = collect(sv.I)
            om = statevec2occmat(hilbert, sv)
            ob = statevec2occbin(hilbert, sv)

            sv1 = occmat2statevec(hilbert, om)
            ob1 = occmat2occbin(hilbert, om)
            @test sv1 == sv
            @test ob1 == ob

            sv2 = occbin2statevec(hilbert, ob)
            om2 = occbin2occmat(hilbert, ob)
            @test sv2 == sv
            @test om2 == om
        end

        @test statevec2locvec(hilbert, [2,1]) == [[1], []]
        @test statevec2locvec(hilbert, [2,2]) == [[1,2], []]
        @test statevec2locvec(hilbert, [2,3]) == [[1,2,2], []]
        @test statevec2locvec(hilbert, [2,4]) == [[1], [2]]
        @test statevec2locvec(hilbert, [2,5]) == [[1,2], [2]]
        @test statevec2locvec(hilbert, [2,6]) == [[1,2,2], [2]]

        @test locvec2statevec(hilbert, [[1], Int[]]) == ([2,1], 1)
        @test locvec2statevec(hilbert, [[1,2], Int[]]) == ([2,2], 1)
        @test locvec2statevec(hilbert, [[1,2,2], Int[]]) == ([2,3], 1)
        @test locvec2statevec(hilbert, [[1], [2]]) == ([2,4], 1)
        @test locvec2statevec(hilbert, [[1,2], [2]]) ==  ([2,5], 1)
        @test locvec2statevec(hilbert, [[1,2,2], [2]]) == ([2,6], 1)

        # TODO(kyungminlee) statevec2occmat, statevec2occbin, locvec2occmat locvec2occbin
    end
end
