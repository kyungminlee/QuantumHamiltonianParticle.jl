using Test
using Particle
using ExactDiagonalization

@testset "particle type" begin
    b = Boson(:b, 3)
    f = Fermion(:f)
    m = HardcoreBoson(:m)

    @testset "constructor" begin
        b1 = Boson("b", 3)
        f1 = Fermion("f")
        m1 = HardcoreBoson("m")
        @test b == b1
        @test f == f1
        @test m == m1
    end

    @test exchangesign(b) == 1
    @test exchangesign(f) == -1
    @test exchangesign(m) == 1

    @test isboson(b)
    @test !isboson(f)
    @test isboson(m)

    @test !isfermion(b)
    @test isfermion(f)
    @test !isfermion(m)

    @test maxoccupancy(b) == 3
    @test maxoccupancy(f) == 1
    @test maxoccupancy(m) == 1

    @test bitwidth(b) == 2
    @test bitwidth(f) == 1
    @test bitwidth(m) == 1
end


@testset "particle sector" begin
    b = Boson(:b, 3)
    f = Fermion(:f)
    m = HardcoreBoson(:m)

    p = ParticleSector(b, f)

    @testset "constructor" begin
        p1 = ParticleSector((b, f))
        @test p == p1
    end

    @test num_particle_species(p) == 2
    @test particle_species(p) == (Boson{:b,3}, Fermion{:f})
    #@show particle_species(p)

    #@test particle_species(p)


end
