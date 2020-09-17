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

        b2 = Boson{:b, 3}()
        f2 = Fermion{:f}()
        m2 = HardcoreBoson{:m}()
        @test b == b2
        @test f == f2
        @test m == m2
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
        p2 = ParticleSector(typeof((b, f)))
        @test p == p2

    end

    @test num_particle_species(p) == 2
    @test particle_species(p) == (Boson{:b,3}, Fermion{:f})

    @test particle_species_name(p, 1) == :b
    @test particle_species_name(p, 2) == :f

    @test exchangesign(p, 1) == 1
    @test exchangesign(p, 2) == -1

    @test bitwidth(p) == 3
    @test bitwidth(p, 1) == 2
    @test bitwidth(p, 2) == 1

    @test bitoffset(p, 1) == 0
    @test bitoffset(p, 2) == 2

    @test bitoffset(p) == [0, 2, 3]


    @test_throws ArgumentError compress(p, [0,0,0], UInt)

    @testset "large" begin
        b1 = Boson(:b, 256)
        p1 = ParticleSector(b1, f)
        @test_throws ArgumentError compress(p1, [0,0], UInt8)
    end

end
