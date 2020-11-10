using Test
using LatticeTools
using ExactDiagonalization
using Particle
using Random


@testset "Particle Projector Operators" begin
    @testset "ParticleProjectorUnitOperator" begin
        p = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0000, 0b0010, 1.0) # should work fine
        @test_throws ArgumentError ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0000, 0b0011, 1.0) # parity bit overlap
        @test_throws ArgumentError ParticleProjectorUnitOperator(0b0101, 0b0110, 0b0000, 0b0000, 1.0) # row
        @test_throws ArgumentError ParticleProjectorUnitOperator(0b0101, 0b0100, 0b1000, 0b0000, 1.0) # row

        zero(p)
    end
end


@testset "make_projector_operator" begin
    p = ParticleSector(Fermion(:f), HardcoreBoson(:b), Spin(:s, 2))
    c(i, j) = ParticleLadderUnit(p, i, j, ANNIHILATION)
    cdag(i, j) = ParticleLadderUnit(p, i, j, CREATION)
    site = ParticleSite([
        ParticleState(p, "__↑", [0,0,0], ( 0,  1)),
        ParticleState(p, "f_↑", [1,0,0], ( 1,  1)),
        ParticleState(p, "f_.", [1,0,1], ( 1,  0)),
        ParticleState(p, "__↓", [0,0,2], ( 0, -1)),
        ParticleState(p, "fb↓", [1,1,2], ( 1, -1)),
    ])
    hilbert_space = ParticleHilbertSpace([site, site, site])
    hsr = represent(hilbert_space)
    rng = MersenneTwister(0)
    for lad in [c(1,2), cdag(1,2), c(2,2), cdag(2,2), c(3,2), cdag(3,2), cdag(1,2)*c(2,1), cdag(1,2)*c(2,1) + cdag(2,2)*0.3]
        opa = embed(hilbert_space, lad)
        opb = make_projector_operator(hilbert_space, lad)
        opc = make_projector_operator(opa)
        for bvec in rand(rng, hsr.basis_list, 5)
            out1a = Dict(collect(get_column_iterator(opa, bvec)))
            out1b = Dict(collect(get_column_iterator(opb, bvec)))
            out1c = Dict(collect(get_column_iterator(opc, bvec)))
            @test out1a == out1b == out1c
            out2a = Dict(collect(get_row_iterator(opa, bvec)))
            out2b = Dict(collect(get_row_iterator(opb, bvec)))
            out2c = Dict(collect(get_row_iterator(opc, bvec)))
            @test out2a == out2b == out2c

            state = SparseState(bvec=>1)
            out3a = SparseState{Float64, UInt}()
            out3b = SparseState{Float64, UInt}()
            out3c = SparseState{Float64, UInt}()
            apply!(out3a, opa, state)
            apply!(out3b, opb, state)
            apply!(out3c, opc, state)
            @test out3a == out3b == out3c
        end
    end

end
