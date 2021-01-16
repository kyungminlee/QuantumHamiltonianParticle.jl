using Test
using Particle
using ExactDiagonalization

@testset "prettyprint" begin

    buf = IOBuffer()
    p = ParticleSector(Fermion(:f), Boson(:b, 3))

    c(i, j) = ParticleLadderUnit(p, i, j, ANNIHILATION)
    cdag(i,j)=ParticleLadderUnit(p, i, j, CREATION)

    @testset "unit" begin
        val = c(1,"xx")
        prettyprint(buf, val)
        @test String(take!(buf)) == "ψ(f,xx)"
        val = cdag(1,"xx")
        prettyprint(buf, val)
        @test String(take!(buf)) == "ψ†(f,xx)"

        val = c(2,"xx")
        prettyprint(buf, val)
        @test String(take!(buf)) == "ψ(b,xx)"

        val = cdag(2,"xx")
        prettyprint(buf, val)
        @test String(take!(buf)) == "ψ†(b,xx)"
    end

    @testset "prod" begin
        val = cdag(1,"x") * c(2,"y")
        prettyprint(buf, val)
        @test String(take!(buf)) == "ψ†(f,x)⋅ψ(b,y)"
    end

    @testset "sum" begin
        val = cdag(1,"x") * c(2,"y")*2 - cdag(1,"z")*1.0
        prettyprint(buf, val)
        @test String(take!(buf)) == "(2.0)⋅ψ†(f,x)⋅ψ(b,y) + (-1.0)⋅ψ†(f,z)"
    end
end
