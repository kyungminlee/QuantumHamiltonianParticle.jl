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

    cdag(i,j) = LadderUnitOperator(p, i, j, CREATION)
    c(i,j) = LadderUnitOperator(p, i, j, ANNIHILATION)

    @testset "product" begin
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
        @test iszero(one(LadderProductOperator{PS, Int, Int})) == false

        @test adjoint(n11) == n11
        @test adjoint(hop) == LadderProductOperator([cdag(1,2), c(1,1)])
        @test ishermitian(n11)
        @test !ishermitian(hop)

        @test exchangesign(c(1,1), c(1,1)) == 1
        @test exchangesign(c(1,1), c(2,1)) == 1
        @test exchangesign(c(2,1), c(2,1)) == -1

        @testset "ordering" begin
            @test c(3,1) < c(1,2) < c(1,1)
            @test !(c(1,1) < c(1,1))
            @test !(c(1,1) < c(1,2))

            @test cdag(1,1) < c(1,1)
            @test cdag(1,1) < c(1,2)

            @test cdag(1,1) < cdag(1,2) < cdag(3,1)

            @test !(cdag(1,1) < cdag(1,1))
            @test !(cdag(1,2) < cdag(1,1))
        end
    end

    @testset "sum" begin
        Σ = LadderSumOperator
        ∏ = LadderProductOperator

        @test Σ(c(1,1)) == Σ([∏([c(1,1)]) => 1])
        n11 = cdag(1,1) * c(1,1)
        n12 = cdag(1,2) * c(1,2)
        hop = cdag(1,2) * c(1,1)

        @test Σ(n11) == Σ([n11=>1])
        @test Σ(n11=>10, n12=>20) == Σ([n11=>10, n12=>20])
        @test Σ(n11=>10, n12=>20) != Σ([n11=>10, n12=>2])
        @test Σ(n11=>10, n12=>20) != Σ([n11=>10])
        @test one(Σ{PS, Int, Int, Float64}) == Σ([one(∏{PS, Int, Int}) => one(Float64)])

        let out = [n11+n12]
            push!(out, n11)
            push!(out, c(1,1))
            @test out[2] == Σ([n11=>1])
            @test out[3] != Σ([n11=>1])
            @test out[3] == Σ([∏([c(1,1)])=>1])
        end

        @test c(1,1)*2.3 == Σ([∏([c(1,1)]) => 2.3])
        @test c(1,1)*2 + c(2,2) == Σ([∏([c(1,1)])=>2, ∏([c(2,2)])=>1])
        @test c(1,1)/2 + c(2,2) == Σ([∏([c(1,1)])=>0.5, ∏([c(2,2)])=>1.0])
        @test c(1,1)//2 + c(2,2) == Σ([∏([c(1,1)])=>1//2, ∏([c(2,2)])=>1//1])
        @test 2*c(1,1) + c(2,2) == Σ([∏([c(1,1)])=>2, ∏([c(2,2)])=>1])
        @test 2\c(1,1) + c(2,2) == Σ([∏([c(1,1)])=>0.5, ∏([c(2,2)])=>1.0])
        @test 2\c(1,1) + c(2,2)*3 == Σ([∏([c(1,1)])=>0.5, ∏([c(2,2)])=>3.0])

        @test n11 + n12 == Σ([n11=>1, n12=>1])
        @test n11 + 2*n12 == Σ([n11=>1, n12=>2])
        @test n11 + n12*2 == Σ([n11=>1, n12=>2])
        @test n11 + n12/2 == Σ([n11=>1.0, n12=>0.5])
        @test n11 + 2\n12 == Σ([n11=>1.0, n12=>0.5])
        @test n11 + n12//2 == Σ([n11=>1//1, n12=>1//2])

        @test adjoint(n11+2*n12) == n11 + 2*n12
        @test ishermitian(2*n12+n11)
        @test !ishermitian(n11 + hop)
        @test ishermitian(hop + adjoint(hop))
    end

    @testset "simplify" begin
        # @show simplify( c(1,1)*c(1,1)*c(1,1) )
        @test simplify( c(2,1)*cdag(2,1)*c(2,1) ) == c(2,1)*1
        @show simplify( c(1,1)*cdag(1,1) + cdag(1,1) * c(1,1) ) == cdag(1,1)*c(1,1) + 1
        @test iszero(simplify( c(2,1)*c(2,1)*c(2,1) ))
    end
end
