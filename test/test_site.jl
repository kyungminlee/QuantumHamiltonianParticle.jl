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

        let state3 = ParticleState(p, "↑f", [0, 1])
            @test state3.name == "↑f"
            @test state3.occupancy_binary == UInt(0b100)
            @test state3.quantum_number == ()
        end

        let state4 = ParticleState(p, "↑f", [0, 1], 1)
            @test state4.name == "↑f"
            @test state4.occupancy_binary == UInt(0b100)
            @test state4.quantum_number == (1,)
        end

        @testset "particle sector" begin
            for s in [state, typeof(state)]
                @test exchangesign(s, 1) == 1
                @test exchangesign(s, 2) == -1
                @test bitwidth(s) == 3
                @test bitwidth(s, 1) == 2
                @test bitwidth(s, 2) == 1
                @test bitoffset(s, 1) == 0
                @test bitoffset(s, 2) == 2
                @test get_bitmask(s, 1) == 0b0011
                @test get_bitmask(s, 2) == 0b0100
                @test compress(s, [0,0]) == 0b000
                @test compress(s, [1,0]) == 0b001
                @test compress(s, [2,0]) == 0b010
                @test compress(s, [0,1]) == 0b100
                @test compress(s, [1,1]) == 0b101
                @test compress(s, [2,1]) == 0b110

                @test extract(s, 0b000) == [0, 0]
                @test extract(s, 0b001) == [1, 0]
                @test extract(s, 0b010) == [2, 0]
                @test extract(s, 0b100) == [0, 1]
                @test extract(s, 0b101) == [1, 1]
                @test extract(s, 0b110) == [2, 1]

                @test numspecies(s) == 2
                @test speciescount(s) == 2
                @test getspecies(s, 1) == Boson{:m, 2}
                @test getspecies(s, 2) == Fermion{:f}
                @test getspeciesname(s, 1) == :m
                @test getspeciesname(s, 2) == :f
            end # testset for
        end # testset particle sector
    end # testset state

    @testset "site" begin
        p = ParticleSector(Boson(:m, 2), Fermion(:f))
        states = [
            ParticleState(p, "↑.", [0, 0], ( 1, 0)),
            ParticleState(p, "0.", [1, 0], ( 0, 0)),
            ParticleState(p, "↑.", [2, 0], (-1, 0)),
            ParticleState(p, "↑f", [0, 1], ( 1, 1)),
            ParticleState(p, "0f", [1, 1], ( 0, 1)),
            ParticleState(p, "↑f", [2, 1], (-1, 1)),
        ]
        site = ParticleSite(states)

        @test qntype(site) == Tuple{Int, Int}
        @test dimension(site) == 6

        for (i, b) in enumerate([0b000, 0b001, 0b010, 0b100, 0b101, 0b110])
            @test get_state_index(site, b) == i
            @test get_state(site, b) == states[i]
        end
        qns = [( 1, 0), ( 0, 0), (-1, 0), ( 1, 1), ( 0, 1), (-1, 1)]
        @test Set(quantum_number_sectors(site)) == Set(qns)
        for (i, qn) in enumerate(qns)
            @test get_quantum_number(site, i) == qn
        end

        @testset "particlesector" begin
            for s in [site, typeof(site)]
                @test exchangesign(s, 1) == 1
                @test exchangesign(s, 2) == -1
                @test bitwidth(s) == 3
                @test bitwidth(s, 1) == 2
                @test bitwidth(s, 2) == 1
                @test bitoffset(s, 1) == 0
                @test bitoffset(s, 2) == 2
                @test get_bitmask(s, 1) == 0b0011
                @test get_bitmask(s, 2) == 0b0100
                @test compress(s, [0,0]) == 0b000
                @test compress(s, [1,0]) == 0b001
                @test compress(s, [2,0]) == 0b010
                @test compress(s, [0,1]) == 0b100
                @test compress(s, [1,1]) == 0b101
                @test compress(s, [2,1]) == 0b110

                @test extract(s, 0b000) == [0, 0]
                @test extract(s, 0b001) == [1, 0]
                @test extract(s, 0b010) == [2, 0]
                @test extract(s, 0b100) == [0, 1]
                @test extract(s, 0b101) == [1, 1]
                @test extract(s, 0b110) == [2, 1]

                @test numspecies(s) == 2
                @test speciescount(s) == 2
                @test getspecies(s, 1) == Boson{:m, 2}
                @test getspecies(s, 2) == Fermion{:f}
                @test getspeciesname(s, 1) == :m
                @test getspeciesname(s, 2) == :f
            end # for s
        end # testset particlesector
    end # testset site
end # state-site
