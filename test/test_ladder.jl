using Test
using LinearAlgebra
using Particle

@testset "ladder" begin
    p = ParticleSector(Boson(:m, 2), Fermion(:f))
    PS = typeof(p)

    @testset "null" begin
        nop1 = ParticleLadderNull(PS)
        nop2 = ParticleLadderNull(p)
        nop3 = ParticleLadderNull{PS}()
        @test nop1 == nop2 == nop3
        @test nop1 === nop2 === nop3
        @test exchangesign(nop1, nop2) == 1
        @test iszero(nop1)
        @test adjoint(nop1) == nop1
        @test ishermitian(nop1)
    end

    @testset "unit" begin
        n = ParticleLadderNull(p)
        b1 = ParticleLadderUnit(p, 1, 2, CREATION)
        @test b1 == ParticleLadderUnit(typeof(p), 1, 2, CREATION)
        @test b1 != ParticleLadderUnit(typeof(p), 1, 2, ANNIHILATION)
        @test b1 != ParticleLadderUnit(typeof(p), 1, 1, CREATION)
        @test b1 != ParticleLadderUnit(typeof(p), 2, 2, CREATION)
        @test !iszero(b1)
        @test adjoint(b1) == ParticleLadderUnit(p, 1, 2, ANNIHILATION)
        @test !ishermitian(b1)
        b2 = ParticleLadderUnit(p, 1, 3, ANNIHILATION)
        c1 = ParticleLadderUnit(p, 2, 1, CREATION)
        c2 = ParticleLadderUnit(p, 2, 2, CREATION)
        @test exchangesign(b1, b2) == 1
        @test exchangesign(b1, c1) == 1
        @test exchangesign(c1, c2) == -1

        @test exchangesign(b1, n) == 1
        @test exchangesign(c1, n) == 1
        @test exchangesign(n, b1) == 1
        @test exchangesign(n, c1) == 1
    end

    cdag(i, j) = ParticleLadderUnit(p, i, j, CREATION)
    c(i, j)    = ParticleLadderUnit(p, i, j, ANNIHILATION)

    @testset "product" begin
        n11 = ParticleLadderProduct([cdag(1,1), c(1,1)])
        n12 = ParticleLadderProduct([cdag(1,2), c(1,2)])
        hop = ParticleLadderProduct([cdag(1,1), c(1,2)])

        @test n11 != n12
        @test n11 == cdag(1,1) * c(1,1)
        @test n11 * c(2,1) == ParticleLadderProduct([cdag(1,1), c(1,1), c(2,1)])
        @test c(2,1) * n11  == ParticleLadderProduct([c(2,1), cdag(1,1), c(1,1)])
        @test n11 * n12 == ParticleLadderProduct([cdag(1,1), c(1,1), cdag(1,2), c(1,2)])

        let
            out = [n11]
            push!(out, c(1,1))
            @test out[2] == ParticleLadderProduct([c(1,1)])
        end

        @test one(ParticleLadderProduct{PS, Int, Int}) == ParticleLadderProduct(ParticleLadderUnit{PS, Int, Int}[])
        @test isone(one(ParticleLadderProduct{PS, Int, Int}))
        @test iszero(one(ParticleLadderProduct{PS, Int, Int})) == false

        @test adjoint(n11) == n11
        @test adjoint(hop) == ParticleLadderProduct([cdag(1,2), c(1,1)])
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
        ∑ = ParticleLadderSum
        ∏ = ParticleLadderProduct

        nop = ParticleLadderNull(p)

        @test ∑(c(1,1)) == ∑([∏([c(1,1)]) => 1])
        n11 = cdag(1,1) * c(1,1)
        n12 = cdag(1,2) * c(1,2)
        hop = cdag(1,2) * c(1,1)

        @test ∑(n11) == ∑([n11=>1])
        @test ∑(n11=>10, n12=>20) == ∑([n11=>10, n12=>20])
        @test ∑(n11=>10, n12=>20) != ∑([n11=>10, n12=>2])
        @test ∑(n11=>10, n12=>20) != ∑([n11=>10])
        @test one(∑{PS, Int, Int, Float64}) == ∑([one(∏{PS, Int, Int}) => one(Float64)])

        let out = [n11+n12]
            push!(out, n11)
            push!(out, c(1,1))
            @test out[2] == ∑([n11=>1])
            @test out[3] != ∑([n11=>1])
            @test out[3] == ∑([∏([c(1,1)])=>1])
        end

        @test c(1,1)*2.3 == ∑([∏([c(1,1)]) => 2.3])

        @test c(1,1) + nop == c(1,1)
        @test nop + c(1,1) == c(1,1)

        @test c(1,1) * nop == nop
        @test nop * c(1,1) == nop

        @test c(1,1) + c(2,2) == ∑([∏([c(1,1)])=>1, ∏([c(2,2)])=>1])
        @test c(1,1)*2 + c(2,2) == ∑([∏([c(1,1)])=>2, ∏([c(2,2)])=>1])
        @test c(1,1)/2 + c(2,2) == ∑([∏([c(1,1)])=>0.5, ∏([c(2,2)])=>1.0])
        @test c(1,1)//2 + c(2,2) == ∑([∏([c(1,1)])=>1//2, ∏([c(2,2)])=>1//1])
        @test 2*c(1,1) + c(2,2) == ∑([∏([c(1,1)])=>2, ∏([c(2,2)])=>1])
        @test 2\c(1,1) + c(2,2) == ∑([∏([c(1,1)])=>0.5, ∏([c(2,2)])=>1.0])
        @test 2\c(1,1) + c(2,2)*3 == ∑([∏([c(1,1)])=>0.5, ∏([c(2,2)])=>3.0])

        @test c(1,1) + n12 == ∑([∏([c(1,1)])=>1, n12=>1])

        @test n11 + n12 == ∑([n11=>1, n12=>1])
        @test n11 + 2*n12 == ∑([n11=>1, n12=>2])
        @test n11 + n12*2 == ∑([n11=>1, n12=>2])
        @test n11 + n12/2 == ∑([n11=>1.0, n12=>0.5])
        @test n11 + 2\n12 == ∑([n11=>1.0, n12=>0.5])
        @test n11 + n12//2 == ∑([n11=>1//1, n12=>1//2])

        @test n11 + c(1,2) == ∑([n11=>1, ∏([c(1,2)])=>1])
        @test c(1,2) + n11 == ∑([∏([c(1,2)])=>1, n11=>1])

        op = n11 + 2*n12

        @test one(ParticleLadderSum{PS, Int, Int, Int}) == one(op)
        @test zero(ParticleLadderSum{PS, Int, Int, Int}) == zero(op)
        @test isone(one(op)) && !iszero(one(op))
        @test iszero(zero(op)) && !isone(zero(op))

        @test adjoint(op) == op
        @test ishermitian(op)
        @test op*2 == 2*n11 + 4*n12
        @test op/2 == 0.5*n11 + 1.0*n12
        @test op//2 == (1//2)*n11 + (1//1)*n12
        @test 2*op == 2*n11 + 4*n12
        @test 2\op == 0.5*n11 + 1.0*n12

        @test c(1,1) * op == ∑([∏([c(1,1), cdag(1,1), c(1,1)])=>1, ∏([c(1,1), cdag(1,2), c(1,2)])=>2])
        @test op * c(1,1) == ∑([∏([cdag(1,1), c(1,1), c(1,1)])=>1, ∏([cdag(1,2), c(1,2), c(1,1)])=>2])

        @test n11 * op == ∑([∏([cdag(1,1), c(1,1), cdag(1,1), c(1,1)])=>1, ∏([cdag(1,1), c(1,1), cdag(1,2), c(1,2)])=>2])
        @test op * n11 == ∑([∏([cdag(1,1), c(1,1), cdag(1,1), c(1,1)])=>1, ∏([cdag(1,2), c(1,2), cdag(1,1), c(1,1)])=>2])

        @test c(1,1) + op == ∑([∏([c(1,1)])=>1, ∏([cdag(1,1), c(1,1)])=>1, ∏([cdag(1,2), c(1,2)])=>2])
        @test op + c(1,1) == ∑([∏([cdag(1,1), c(1,1)])=>1, ∏([cdag(1,2), c(1,2)])=>2, ∏([c(1,1)])=>1, ])

        @test n11 + op == ∑([n11=>1, ∏([cdag(1,1), c(1,1)])=>1, ∏([cdag(1,2), c(1,2)])=>2])
        @test op + n11 == ∑([∏([cdag(1,1), c(1,1)])=>1, ∏([cdag(1,2), c(1,2)])=>2, n11=>1])

        @test c(1,1) + 1 == ∑([∏([c(1,1)])=>1, one(∏{PS, Int, Int})=>1])
        @test 1 + c(1,1) == ∑([one(∏{PS, Int, Int})=>1, ∏([c(1,1)])=>1])

        @test c(1,1) - 1 == ∑([∏([c(1,1)])=>1, one(∏{PS, Int, Int})=>-1])
        @test 1 - c(1,1) == ∑([one(∏{PS, Int, Int})=>1, ∏([c(1,1)])=>-1])

        @test n11 + 1 == ∑([n11=>1, one(n11)=>1])
        @test 1 + n11 == ∑([one(n11)=>1, n11=>1])

        @test op + 1 == ∑([∏([cdag(1,1), c(1,1)])=>1, ∏([cdag(1,2), c(1,2)])=>2, one(∏{PS, Int, Int})=>1 ])
        @test 1 + op == ∑([ one(∏{PS, Int, Int})=>1, ∏([cdag(1,1), c(1,1)])=>1, ∏([cdag(1,2), c(1,2)])=>2 ])

        @test !ishermitian(n11 + hop)
        @test ishermitian(hop + adjoint(hop))
    end

    @testset "unary operation" begin
        ∑ = ParticleLadderSum
        ∏ = ParticleLadderProduct

        n = ParticleLadderNull(p)

        @test n == -n
        @test +c(1,1) == c(1,1)
        @test -c(1,1) == ∑([∏([c(1,1)])=>-1])
        n11 = cdag(1,1)*c(1,1)
        @test +n11 == n11
        @test -n11 == ∑([n11=>-1])
        hop = cdag(1,2)*c(1,1) + cdag(1,1)*c(1,2)
        @test +hop == hop
        @test -hop == ∑([cdag(1,2)*c(1,1)=>-1, cdag(1,1)*c(1,2)=>-1])
    end


    @testset "simplify" begin
        @test simplify( c(2,1) ) == c(2,1)
        @test simplify( c(2,1)*cdag(2,1)*c(2,1) ) == c(2,1)*1
        @test iszero(simplify( c(2,1)*cdag(2,1) + cdag(2,1)*c(2,1) + (-1)))
        @test iszero(simplify( c(2,1)*cdag(2,1) + cdag(2,1)*c(2,1) - 1))
        @test iszero(simplify( c(2,1)*c(2,1)*c(2,1) ))
    end

    @testset "with Hilbert" begin
        # p = ParticleSector(Boson(:m, 2), Fermion(:f))

        site = ParticleSite([
            ParticleState(p, "__", [0, 0], (0, 0)),
            ParticleState(p, "b_", [1, 0], (1, 0)),
            ParticleState(p, "B_", [2, 0], (2, 0)),
            ParticleState(p, "_f", [0, 1], (0, 1)),
            ParticleState(p, "bf", [1, 1], (1, 1)),
            ParticleState(p, "Bf", [2, 1], (2, 1)),
        ])
        hs = ParticleHilbertSpace([site, site, site])
        @test get_fermion_parity(hs, cdag(1, 3), 0b000_000_001) == 0
        @test get_fermion_parity(hs, cdag(2, 3), 0b000_000_001) == 0

        @test get_fermion_parity(hs, cdag(1, 3), 0b000_000_100) == 0
        @test get_fermion_parity(hs, cdag(2, 3), 0b000_000_100) == 1

        @test get_fermion_parity(hs, cdag(1, 1), 0b001_000_000) == 0
        @test get_fermion_parity(hs, cdag(2, 1), 0b001_000_000) == 0

        @test get_fermion_parity(hs, cdag(1, 1), 0b100_000_000) == 0
        @test get_fermion_parity(hs, cdag(2, 1), 0b100_000_000) == 0
    end
end
