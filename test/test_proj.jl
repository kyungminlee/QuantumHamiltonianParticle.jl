using Test
using LatticeTools
using ExactDiagonalization
using Particle
using Random

@testset "ParticleProjectorUnitOperator" begin
    p = ParticleSector(Fermion(:f), HardcoreBoson(:b), Spin(:s, 1))
    c(i, j) = ParticleLadderUnit(p, i, j, ANNIHILATION)
    cdag(i, j) = ParticleLadderUnit(p, i, j, CREATION)
    site = ParticleSite([
        ParticleState(p, "__↑", [0,0,0], ( 0,  1)),
        ParticleState(p, "f_↑", [1,0,0], ( 1,  1)),
        ParticleState(p, "__↓", [0,0,1], ( 0, -1)),
        ParticleState(p, "fb↓", [1,1,1], ( 1, -1)),
    ])
    hilbert_space = ParticleHilbertSpace([site, site, site])
    hsr = represent(hilbert_space)
    rng = MersenneTwister(0)
    for lad in [c(1,2), cdag(1,2), c(2,2), cdag(2,2), c(3,2), cdag(3,2), cdag(1,2)*c(2,1)]
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
