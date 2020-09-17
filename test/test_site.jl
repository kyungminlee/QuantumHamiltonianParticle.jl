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

        let
            state3 = ParticleState(p, "↑f", [0, 1])
            @test state3.name == "↑f"
            @test state3.occupancy_binary == UInt(0b100)
            @test state3.quantum_number == ()
        end

        let
            state4 = ParticleState(p, "↑f", [0, 1], 1)
            @test state4.name == "↑f"
            @test state4.occupancy_binary == UInt(0b100)
            @test state4.quantum_number == (1,)
        end

        @test exchangesign(state, 1) == 1
        @test exchangesign(state, 2) == -1
        @test bitwidth(state) == 3
        @test bitwidth(state, 1) == 2
        @test bitwidth(state, 2) == 1
        @test bitoffset(state, 1) == 0
        @test bitoffset(state, 2) == 2
        @test get_bitmask(state, 1) == 0b0011
        @test get_bitmask(state, 2) == 0b0100
        @test compress(state, [0,0]) == 0b000
        @test compress(state, [1,0]) == 0b001
        @test compress(state, [2,0]) == 0b010
        @test compress(state, [0,1]) == 0b100
        @test compress(state, [1,1]) == 0b101
        @test compress(state, [2,1]) == 0b110

        @test extract(state, 0b000) == [0, 0]
        @test extract(state, 0b001) == [1, 0]
        @test extract(state, 0b010) == [2, 0]
        @test extract(state, 0b100) == [0, 1]
        @test extract(state, 0b101) == [1, 1]
        @test extract(state, 0b110) == [2, 1]

        @test numspecies(state) == 2
        @test speciescount(state) == 2
        @test getspecies(state, 1) == Boson{:m, 2}
        @test getspecies(state, 2) == Fermion{:f}
        @test getspeciesname(state, 1) == :m
        @test getspeciesname(state, 2) == :f
    end

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




        @test exchangesign(site, 1) == 1
        @test exchangesign(site, 2) == -1
        @test bitwidth(site) == 3
        @test bitwidth(site, 1) == 2
        @test bitwidth(site, 2) == 1
        @test bitoffset(site, 1) == 0
        @test bitoffset(site, 2) == 2
        @test get_bitmask(site, 1) == 0b0011
        @test get_bitmask(site, 2) == 0b0100
        @test compress(site, [0,0]) == 0b000
        @test compress(site, [1,0]) == 0b001
        @test compress(site, [2,0]) == 0b010
        @test compress(site, [0,1]) == 0b100
        @test compress(site, [1,1]) == 0b101
        @test compress(site, [2,1]) == 0b110

        @test extract(site, 0b000) == [0, 0]
        @test extract(site, 0b001) == [1, 0]
        @test extract(site, 0b010) == [2, 0]
        @test extract(site, 0b100) == [0, 1]
        @test extract(site, 0b101) == [1, 1]
        @test extract(site, 0b110) == [2, 1]

        @test numspecies(site) == 2
        @test speciescount(site) == 2
        @test getspecies(site, 1) == Boson{:m, 2}
        @test getspecies(site, 2) == Fermion{:f}
        @test getspeciesname(site, 1) == :m
        @test getspeciesname(site, 2) == :f


    end
end
