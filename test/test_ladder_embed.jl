using Test
using LinearAlgebra
using Random

using ExactDiagonalization
using LatticeTools
using Particle

@testset "ParticleLadderEmbedding" begin
    p = ParticleSector(Boson(:m, 2), Fermion(:f))
    cdag(i, j) = ParticleLadderUnit(p, i, j, CREATION)
    c(i, j)    = ParticleLadderUnit(p, i, j, ANNIHILATION)

    site1 = ParticleSite([
        ParticleState(p, "__", [0, 0], (0, 0)),
        ParticleState(p, "b_", [1, 0], (1, 0)),
        ParticleState(p, "B_", [2, 0], (2, 0)),
        ParticleState(p, "_f", [0, 1], (0, 1)),
        ParticleState(p, "bf", [1, 1], (1, 1)),
        ParticleState(p, "Bf", [2, 1], (2, 1)),
    ])
    hs1 = ParticleHilbertSpace([site1, site1, site1])

    @testset "Constructor & Comparison" begin
        site2 = ParticleSite([
            ParticleState(p, "__", [0, 0], (0, 0)),
            ParticleState(p, "b_", [1, 0], (1, 0)),
            ParticleState(p, "_f", [0, 1], (0, 1)),
            ParticleState(p, "bf", [1, 1], (1, 1)),
        ])
        hs2 = ParticleHilbertSpace([site2, site2])
        @test cdag(1,1) == cdag(1,1)
        @test embed(hs1, cdag(1,1)) == embed(hs1, cdag(1,1))
        @test embed(hs1, cdag(1,1)) != embed(hs1, cdag(1,2))
        @test embed(hs1, cdag(1,1)) != embed(hs2, cdag(1,1))

        @test embed(hs1, cdag(1,1) + cdag(1,2)) != embed(hs1, cdag(1,2) + cdag(1,1))
        @test isequiv(embed(hs1, cdag(1,1) + cdag(1,2)), embed(hs1, cdag(1,2) + cdag(1,1)))

        # TODO(2020-11-09): implement isapprox
        # @test embed(hs1, cdag(1,1) + (1+1E-12) * cdag(1,2)) != embed(hs1, cdag(1,1) + cdag(1,2))
        # @test isapprox(embed(hs1, cdag(1,1) + (1+1E-12) * cdag(1,2)), embed(hs1, cdag(1,1) + 1.0 * cdag(1,2)); atol=1E-8, rtol=1E-8)
    end

    @testset "Methods" begin
    end

    @testset "Ladder Iterator" begin
        p = ParticleSector(Boson(:m, 2), Fermion(:f))
        cdag(i, j) = ParticleLadderUnit(p, i, j, CREATION)
        c(i, j)    = ParticleLadderUnit(p, i, j, ANNIHILATION)
        site = ParticleSite([
            ParticleState(p, "__", [0, 0], (0, 0)),
            ParticleState(p, "b_", [1, 0], (1, 0)),
            ParticleState(p, "B_", [2, 0], (2, 0)),
            ParticleState(p, "_f", [0, 1], (0, 1)),
            ParticleState(p, "bf", [1, 1], (1, 1)),
            ParticleState(p, "Bf", [2, 1], (2, 1)),
        ])
        hs = ParticleHilbertSpace([site, site, site])

        # ParticleLadderNull
        nop = ParticleLadderNull(p)
        @test isempty(collect(get_column_iterator(hs, nop, UInt(0b000_000_000))))
        @test isempty(collect(get_row_iterator(hs, nop, UInt(0b000_000_000))))

        hsr = represent(hs)
        rng = MersenneTwister(0)
        for op in [cdag(1,2), cdag(2,2), c(1,2), c(2,2), cdag(1,3)*c(1,1), cdag(2,3) * c(2,1), cdag(2,3)*c(2,1)*2 + c(1,1)*0.5]
            for bvec in rand(rng, hsr.basis_list, 8)
                @test collect(get_column_iterator(hs, op, bvec)) == collect(get_column_iterator(embed(hs, op), bvec))
                @test collect(get_row_iterator(hs, op, bvec)) == collect(get_row_iterator(embed(hs, op), bvec))
            end
            for brow in rand(rng, hsr.basis_list, 3), bcol in rand(rng, hsr.basis_list, 3)
                @test get_element(hs, op, brow, bcol) == get_element(embed(hs, op), brow, bcol)
            end
        end
    end

    @testset "Unary Operation" begin
        op = embed(hs1, cdag(1,2))
        @test +op == op
        @test -op == embed(hs1, -cdag(1,2))
        @test adjoint(op) == embed(hs1, c(1,2))
        @test !ishermitian(op)
        @test ishermitian(embed(hs1, cdag(1,2)*c(1,2)))
        @test !ishermitian(embed(hs1, cdag(1,2)*c(1,1)))
        @test ishermitian(embed(hs1, cdag(1,2)*c(1,1) + cdag(1,1)*c(1,2)))
    end

    @testset "Binary Operation" begin
        site2 = ParticleSite([
            ParticleState(p, "__", [0, 0], (0, 0)),
            ParticleState(p, "b_", [1, 0], (1, 0)),
            ParticleState(p, "_f", [0, 1], (0, 1)),
            ParticleState(p, "bf", [1, 1], (1, 1)),
        ])
        hs2 = ParticleHilbertSpace([site2, site2])

        op1 = embed(hs1, cdag(1,2))
        op2 = embed(hs2, c(1,1))
        op3 = embed(hs1, c(1,1))
        @test_throws ArgumentError op1 + op2
        @test_throws ArgumentError op1 - op2
        @test_throws ArgumentError op1 * op2

        @test op1 + op3 == embed(hs1, cdag(1,2) + c(1,1))
        @test op1 - op3 == embed(hs1, cdag(1,2) - c(1,1))
        @test op1 * op3 == embed(hs1, cdag(1,2) * c(1,1))

        @test op1 * 2.0 == embed(hs1, cdag(1,2) * 2.0)
        @test 2.0 * op1 == embed(hs1, 2.0 * cdag(1,2))
        @test op1 / 2.0 == embed(hs1, cdag(1,2) / 2.0)
        @test op1 // 2 == embed(hs1, cdag(1,2) // 2)
        @test 2.0 \ op1 == embed(hs1, 2.0 \ cdag(1,2))
    end

    @testset "apply!" begin
    end

    @testset "simplify" begin
    end
end
