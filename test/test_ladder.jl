using Test
using LinearAlgebra
using Particle

@testset "ladder" begin
    p = ParticleSector(Boson(:m, 2), Fermion(:f))
    PS = typeof(p)
    @testset "unit" begin
        c = LadderUnitOperator(p, 1, 2, CREATION)
        @test c == LadderUnitOperator(typeof(p), 1, 2, CREATION)
        @test c != LadderUnitOperator(typeof(p), 1, 2, ANNIHILATION)
        @test c != LadderUnitOperator(typeof(p), 1, 1, CREATION)
        @test c != LadderUnitOperator(typeof(p), 2, 2, CREATION)
        @test !iszero(c)
        @test adjoint(c) == LadderUnitOperator(p, 1, 2, ANNIHILATION)
        @test !ishermitian(c)
    end

    @testset "product" begin
        cdag(i,j) = LadderUnitOperator(p, i, j, CREATION)
        c(i,j) = LadderUnitOperator(p, i, j, ANNIHILATION)
        n11 = LadderProductOperator([cdag(1,1), c(1,1)])
        n12 = LadderProductOperator([cdag(1,2), c(1,2)])
        hop = LadderProductOperator([cdag(1,1), c(1,2)])

        @test n11 != n12
        @test n11 == cdag(1,1) * c(1,1)
        @test n11 * c(2,1) == LadderProductOperator([cdag(1,1), c(1,1), c(2,1)])
        @test c(2,1) * n11  == LadderProductOperator([c(2,1), cdag(1,1), c(1,1)])
        @test n11 * n12 == LadderProductOperator([cdag(1,1), c(1,1), cdag(1,2), c(1,2)])

        let
            out = [n11]
            push!(out, c(1,1))
            @test out[2] == LadderProductOperator([c(1,1)])
        end

        @test one(LadderProductOperator{PS, Int, Int}) == LadderProductOperator(LadderUnitOperator{PS, Int, Int}[])
        @test isone(one(LadderProductOperator{PS, Int, Int}))
        @test adjoint(n11) == n11
        @test adjoint(hop) == LadderProductOperator([cdag(1,2), c(1,1)])
        @test ishermitian(n11)
        @test !ishermitian(hop)

    end
end
