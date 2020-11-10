using Test
using LatticeTools
using ExactDiagonalization
using Particle
using Random


@testset "Particle Projector Operators" begin
    @testset "ParticleProjectorUnitOperator" begin
        pi = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1) # should work fine
        pd = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1.0) # should work fine
        pz = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1.0 + 0.0im) # should work fine
        @test_throws ArgumentError ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0000, 0b0011, 1.0) # parity bit overlap
        @test_throws ArgumentError ParticleProjectorUnitOperator(0b0101, 0b0110, 0b0000, 0b0000, 1.0) # row
        @test_throws ArgumentError ParticleProjectorUnitOperator(0b0101, 0b0100, 0b1000, 0b0000, 1.0) # row

        @test pi == ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1)
        @test pi != ParticleProjectorUnitOperator(0b1101, 0b0100, 0b0001, 0b0010, 1)
        @test pi != ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0001, 0b0010, 1)
        @test pi != ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0000, 0b0010, 1)
        @test pi != ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0000, 1)
        @test pi != ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 2)
        @test isapprox(pi, ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1))
        @test isapprox(pd, ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1))
        @test isapprox(pd, ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1.0))
        @test isapprox(pd, ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1.0 + 1E-12))
        @test isapprox(pd, ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1.0 + 1E-12im))

        x = zero(pi)
        @test x.bitmask == 0x0
        @test x.bitrow == 0x0
        @test x.bitcol == 0x0
        @test x.parity_bitmask == 0x0
        @test iszero(x.amplitude)
        @test typeof(x.amplitude) == Int

        x = zero(pd)
        @test x.bitmask == 0x0
        @test x.bitrow == 0x0
        @test x.bitcol == 0x0
        @test x.parity_bitmask == 0x0
        @test iszero(x.amplitude)
        @test typeof(x.amplitude) == Float64

        x = zero(pz)
        @test x.bitmask == 0x0
        @test x.bitrow == 0x0
        @test x.bitcol == 0x0
        @test x.parity_bitmask == 0x0
        @test iszero(x.amplitude)
        @test typeof(x.amplitude) == ComplexF64

        x = one(pi)
        @test x.bitmask == 0x0
        @test x.bitrow == 0x0
        @test x.bitcol == 0x0
        @test x.parity_bitmask == 0x0
        @test isone(x.amplitude)
        @test typeof(x.amplitude) == Int

        x = one(pd)
        @test x.bitmask == 0x0
        @test x.bitrow == 0x0
        @test x.bitcol == 0x0
        @test x.parity_bitmask == 0x0
        @test isone(x.amplitude)
        @test typeof(x.amplitude) == Float64

        x = one(pz)
        @test x.bitmask == 0x0
        @test x.bitrow == 0x0
        @test x.bitcol == 0x0
        @test x.parity_bitmask == 0x0
        @test isone(x.amplitude)
        @test typeof(x.amplitude) == ComplexF64
    end

    @testset "multiplication" begin
        @testset "mismatch" begin
            p1 = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 2)
            p2 = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 3)
            @test iszero(p1 * p2)
        end
        @testset "match" begin
            # fermion parity noflip
            p1 = ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3)
            p2 = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0000, 2)
            @test p1 * p2 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0001, 0b0010, 6)

            # fermion parity reflip
            p3 = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 2)
            @test p1 * p3 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0001, 0b0000, 6)
        end
    end

    @testset "scale" begin
        p1 = ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3)
        @test p1 * 4 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 12)
        @test 4 * p1 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 12)
        @test p1 / 4 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 0.75)
        @test p1 // 4 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3//4)
        @test 4 \ p1 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 0.75)
        @test -p1 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, -3)
        @test +p1 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3)
        n = 2
        @test p1^n == p1 * p1
        @test iszero(p1^n)

        p2 = ParticleProjectorUnitOperator(0b0101, 0b0101, 0b0101, 0b0010, 3)
        n = 6
        @test p2^n == p2 * p2 * p2 * p2 * p2 * p2
        @test !iszero(p2^6)
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
