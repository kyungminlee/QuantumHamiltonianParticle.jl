using Test
using Random
using LinearAlgebra

using LatticeTools
using ExactDiagonalization
using Particle


@testset "Particle Projector Operators" begin
    @testset "ParticleProjectorUnitOperator" begin
        @testset "constructor" begin
            p0 = ParticleProjectorUnitOperator{UInt, Float64}(10)
            @test isa(p0, ParticleProjectorUnitOperator{UInt, Float64})
            @test p0.bitmask == 0x0
            @test p0.bitrow == 0x0
            @test p0.bitcol == 0x0
            @test p0.parity_bitmask == 0x0
            @test p0.amplitude == 10.0

            pi = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1) # should work fine
            pd = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1.0) # should work fine
            pz = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1.0 + 0.0im) # should work fine
            @test_throws ArgumentError ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0000, 0b0011, 1.0) # parity bit overlap
            @test_throws ArgumentError ParticleProjectorUnitOperator(0b0101, 0b0110, 0b0000, 0b0000, 1.0) # row
            @test_throws ArgumentError ParticleProjectorUnitOperator(0b0101, 0b0100, 0b1000, 0b0000, 1.0) # row

            ParticleProjectorUnitOperator{UInt, Float64}(0b0101, 0b0100, 0b0001, 0b0010, 1) # should work fine
            ParticleProjectorUnitOperator{UInt, Float64}(0b0101, 0b0100, 0b0001, 0b0010, 1.0) # should work fine
            ParticleProjectorUnitOperator{UInt, ComplexF64}(0b0101, 0b0100, 0b0001, 0b0010, 1.0 + 0.0im) # should work fine
            @test_throws ArgumentError ParticleProjectorUnitOperator{UInt, Float64}(0b0101, 0b0100, 0b0000, 0b0011, 1.0) # parity bit overlap
            @test_throws ArgumentError ParticleProjectorUnitOperator{UInt, Float64}(0b0101, 0b0110, 0b0000, 0b0000, 1.0) # row
            @test_throws ArgumentError ParticleProjectorUnitOperator{UInt, Float64}(0b0101, 0b0100, 0b1000, 0b0000, 1.0) # row

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

        @testset "scale and unary" begin
            p1 = ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3)
            @test p1 * 4 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 12)
            @test 4 * p1 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 12)
            @test p1 / 4 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 0.75)
            @test p1 // 4 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3//4)
            @test 4 \ p1 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 0.75)
            @test -p1 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, -3)
            @test +p1 == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3)

            p2 = ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3.0 + 4.0im)
            @test isa(real(p2), ParticleProjectorUnitOperator{UInt8, Float64})
            @test real(p2) == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3.0)
            @test isa(imag(p2), ParticleProjectorUnitOperator{UInt8, Float64})
            @test imag(p2) == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 4.0)
            @test isa(conj(p2), ParticleProjectorUnitOperator{UInt8, ComplexF64})
            @test conj(p2) == ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3.0-4.0im)
            @test isa(adjoint(p2), ParticleProjectorUnitOperator{UInt8, ComplexF64})
            @test adjoint(p2) == ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0000, 0b0010, 3.0-4.0im)
            @test isa(transpose(p2), ParticleProjectorUnitOperator{UInt8, ComplexF64})
            @test transpose(p2) == ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0000, 0b0010, 3.0+4.0im)
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

        @testset "power" begin
            p1 = ParticleProjectorUnitOperator(0b0101, 0b0000, 0b0100, 0b0010, 3)
            n = 2
            @test p1^n == p1 * p1
            @test iszero(p1^n)

            p2 = ParticleProjectorUnitOperator(0b0101, 0b0101, 0b0101, 0b0010, 3)
            n = 6
            @test p2^n == p2 * p2 * p2 * p2 * p2 * p2
            @test !iszero(p2^6)
            @test (p2^6).amplitude == 3^6
        end
    end

    @testset "ParticleProjectorSumOperator" begin
        P∑ = ParticleProjectorSumOperator
        PU = ParticleProjectorUnitOperator

        @testset "zero and one" begin
            z = zero(P∑{UInt, Float64})
            @test isa(z, P∑{UInt, Float64})
            @test z == P∑(PU{UInt, Float64}[])
            o = one(P∑{UInt, Float64})
            @test isa(o, P∑{UInt, Float64})
            @test o == P∑([one(PU{UInt, Float64})])
        end

        @testset "equality" begin
            p1 = PU(0b0101, 0b0100, 0b0001, 0b0010, 2)
            p2 = PU(0b0101, 0b0001, 0b0100, 0b0010, 2.0+0.0im)
            p3 = PU(0b0101, 0b0001, 0b0100, 0b0010, 2.0-0.0im)
            @test P∑([p1, p2]) == P∑([p1, p3])
            @test P∑([p1, p2]) != P∑([p2, p1])  # Not necessary though.
        end

        @testset "unary functions" begin
            p1 = PU(0b0101, 0b0100, 0b0001, 0b0010, 2.0+im)
            p2 = PU(0b0101, 0b0001, 0b0100, 0b0010, 2.0+im) # no conjugation
            p3 = PU(0b0101, 0b0001, 0b0100, 0b0010, 2.0-im) # conjugation
            s1 = P∑([p1, p2])
            @test !ishermitian(s1)
            s2 = P∑([p1, p3])
            @test ishermitian(s2)

            @test +s1 == s1
            @test -s1 == P∑([-p1, -p2])

            @test isa(real(s2), P∑{UInt8, Float64})
            @test isa(imag(s2), P∑{UInt8, Float64})
            @test isa(conj(s2), P∑{UInt8, ComplexF64})
            @test isa(transpose(s2), P∑{UInt8, ComplexF64})
            @test isa(adjoint(s2), P∑{UInt8, ComplexF64})

            @test real(s2) == P∑([PU(0b0101, 0b0100, 0b0001, 0b0010, 2.0), PU(0b0101, 0b0001, 0b0100, 0b0010, 2.0)])
            @test imag(s2) == P∑([PU(0b0101, 0b0100, 0b0001, 0b0010, 1.0), PU(0b0101, 0b0001, 0b0100, 0b0010,-1.0)])
            @test conj(s2) == P∑([PU(0b0101, 0b0100, 0b0001, 0b0010, 2.0-1.0im), PU(0b0101, 0b0001, 0b0100, 0b0010, 2.0+1.0im)])
            @test adjoint(s2) == P∑([PU(0b0101, 0b0001, 0b0100, 0b0010, 2.0-1.0im), PU(0b0101, 0b0100, 0b0001, 0b0010, 2.0+1.0im)])
            @test transpose(s2) == P∑([PU(0b0101, 0b0001, 0b0100, 0b0010, 2.0+1.0im), PU(0b0101, 0b0100, 0b0001, 0b0010, 2.0-1.0im)])
        end

        @testset "binary functions" begin
            p0 = PU(0b0000, 0b0000, 0b0000, 0b0000, 1)
            p1 = PU(0b0101, 0b0100, 0b0100, 0b0010, 2)
            p2 = PU(0b0101, 0b0001, 0b0100, 0b0010, 3 + 4im)
            p3 = PU(0b0101, 0b0100, 0b0001, 0b0010, 5.0)

            s1 = P∑([p1, p2])
            s2 = P∑([p1, p3])

            @test s1 * 2 == P∑([p1*2, p2*2])
            @test 2 * s1 == P∑([2*p1, 2*p2])
            @test s1 * (2.5+0.5im) == P∑([p1*(2.5+0.5im), p2*(2.5+0.5im)])
            @test (2.5+0.5im) * s1 == P∑([(2.5+0.5im)*p1, (2.5+0.5im)*p2])

            @test s1 / 2 == P∑([p1 / 2, p2 / 2])
            @test s1 // 2 == P∑([p1 // 2, p2 // 2])
            @test 2 \ s1 == P∑([2\p1, 2\p2])
            @test s1 / (2.5+0.5im) == P∑([p1 / (2.5+0.5im), p2 / (2.5+0.5im)])
            @test (2.5+0.5im) \ s1 == P∑([(2.5+0.5im) \ p1, (2.5+0.5im) \ p2])

            @test 10 + s1 == P∑([p0*10, p1, p2])
            @test s1 + 10.0 == P∑([p1, p2, p0*10.0])

            @test p1 + s1 == P∑([p1, p1, p2])
            @test p1 - s1 == P∑([p1, -p1, -p2])

            @test s1 + p1 == P∑([p1, p2, p1])
            @test s1 - p1 == P∑([p1, p2, -p1])

            @test s1 + s2 == P∑([p1, p2, p1, p3])
            @test s1 - s2 == P∑([p1, p2, -p1, -p3])

            @test s1 * s2 == P∑([p1*p1, p1*p3, p2*p1, p2*p3]) # Not necessarily true
        end
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


@testset "simplify" begin
    p0 = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 0)
    p1 = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1)
    p2 = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0000, 10)
    p3 = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 100.0)
    p4 = ParticleProjectorUnitOperator(0b0000, 0b0000, 0b0000, 0b0000, 1000.0)

    q0 = simplify(p0)
    @test isa(q0, NullOperator)
    q1 = ParticleProjectorSumOperator([p1, p2, p3, p4])
    q2 = simplify(q1)

    @test q1 != ParticleProjectorSumOperator([
        ParticleProjectorUnitOperator(0b0000, 0b0000, 0b0000, 0b0000, 1000.0),
        ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0000, 10.0),
        ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 101.0),
    ])
    @test q2 == ParticleProjectorSumOperator([
        ParticleProjectorUnitOperator(0b0000, 0b0000, 0b0000, 0b0000, 1000.0),
        ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0000, 10.0),
        ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 101.0),
    ])
    q3 = simplify(ParticleProjectorSumOperator([p1, -p1]))
    @test iszero(q3)
    @test q3 == NullOperator()

    q4 = simplify(ParticleProjectorSumOperator([p1, p1, -p1]))
    @test typeof(q4) == ParticleProjectorUnitOperator{UInt8, Int}
    @test q4 == p1

    q5 = simplify(ParticleProjectorSumOperator([p1, p1, p4, -p1]))
    @test typeof(q5) == ParticleProjectorSumOperator{UInt8, Float64}
    @test q5 == ParticleProjectorSumOperator([p4, p1])

    p6 = ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1.0 + 1E-12im)
    q6 = simplify(p6)
    @test q6 == ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1.0)
    @test typeof(q6) == ParticleProjectorUnitOperator{UInt8, Float64}

    p7 = ParticleProjectorSumOperator([p6])
    q7 = simplify(p7)
    @test q7 == ParticleProjectorUnitOperator(0b0101, 0b0100, 0b0001, 0b0010, 1.0)
    @test typeof(q7) == ParticleProjectorUnitOperator{UInt8, Float64}

end
