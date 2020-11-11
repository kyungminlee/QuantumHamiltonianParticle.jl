using Test
using Particle

@testset "Promotion" begin
    p = ParticleSector(Fermion(:f), HardcoreBoson(:b))
    PS = typeof(p)
    p2 = ParticleSector(Fermion(:f), Boson(:m, 2))
    nop = ParticleLadderNull(p)
    nop2 = ParticleLadderNull(p2)

    c(i, j) = ParticleLadderUnit(p, i, j, ANNIHILATION)
    cdag(i, j) = ParticleLadderUnit(p, i, j, CREATION)

    c2(i, j) = ParticleLadderUnit(p2, i, j, ANNIHILATION)

    @testset "Product" begin
        out = [c(1,2), c(1,1)*cdag(1,2)]
        push!(out, c(1,1))
        push!(out, c(1,1)*c(2,1))
        @test_throws MethodError push!(out, c2(1,1))
        @test_throws MethodError push!(out, c2(1,1)*c2(1,1))
    end

    @testset "Sum" begin
        out = [nop, c(1,2), c(1,1)*cdag(1,2), c(1,2)+cdag(1,2), 1.2 * c(1,1)]
        push!(out, c(1,1))
        push!(out, c(1,1)*c(2,1))
        push!(out, c(1,1)*c(2,1) + c(1,1))
        @test_throws MethodError push!(out, c2(1,1))
        @test_throws MethodError push!(out, c2(1,1)*c2(1,1))

        out = [c(1,2)+cdag(1,2)]
        @test_throws InexactError push!(out, 1.2 * c(1,1))
    end
end
