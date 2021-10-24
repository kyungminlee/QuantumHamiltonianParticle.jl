using Test
using QuantumHamiltonianParticle
using QuantumHamiltonian

@testset "particle type" begin
    b = Boson(:b, 3)
    f = Fermion(:f)
    m = HardcoreBoson(:m)
    s = Spin("S=1", 3)

    @testset "constructor" begin
        b1 = Boson("b", 3)
        f1 = Fermion("f")
        m1 = HardcoreBoson("m")
        s1 = Spin("S=1", 3)
        s1p = Spin(Symbol("S=1"), 3)
        @test b == b1
        @test f == f1
        @test m == m1
        @test s == s1
        @test s == s1p

        b2 = Boson{:b, 3}()
        f2 = Fermion{:f}()
        m2 = HardcoreBoson{:m}()
        s2 = Spin{Symbol("S=1"), 3}()
        @test b == b2
        @test f == f2
        @test m == m2
        @test s == s2
    end

    @test exchangesign(b) == 1
    @test exchangesign(f) == -1
    @test exchangesign(m) == 1
    @test exchangesign(s) == 1

    @test isboson(b)
    @test !isboson(f)
    @test isboson(m)
    @test !isboson(s)

    @test !isfermion(b)
    @test isfermion(f)
    @test !isfermion(m)
    @test !isfermion(s)

    @test !isspin(b)
    @test !isspin(f)
    @test !isspin(m)
    @test isspin(s)

    @test maxoccupancy(b) == 3
    @test maxoccupancy(f) == 1
    @test maxoccupancy(m) == 1
    @test maxoccupancy(s) == 3

    @test bitwidth(b) == 2
    @test bitwidth(f) == 1
    @test bitwidth(m) == 1
    @test bitwidth(s) == 2
end


@testset "particle sector" begin
    b = Boson(:b, 5)
    f = Fermion(:f)
    m = HardcoreBoson(:m)
    s = Spin("S=3/2", 4)

    p = ParticleSector(b, f)

    @testset "constructor" begin
        p1 = ParticleSector((b, f))
        @test p == p1
        p2 = ParticleSector(typeof((b, f)))
        @test p == p2
    end

    @test numspecies(p) == 2
    @test speciescount(p) == 2
    @test getspecies(p) == (Boson{:b,5}, Fermion{:f})

    @test getspeciesname(p, 1) == :b
    @test getspeciesname(p, 2) == :f

    @test exchangesign(p, 1) == 1
    @test exchangesign(p, 2) == -1

    @test exchangesign(p, 1, 1) == 1
    @test exchangesign(p, 1, 2) == 1
    @test exchangesign(p, 2, 1) == 1
    @test exchangesign(p, 2, 2) == -1

    @test bitwidth(p) == 4
    @test bitwidth(p, 1) == 3
    @test bitwidth(p, 2) == 1

    @test bitoffset(p, 1) == 0
    @test bitoffset(p, 2) == 3

    @test bitoffset(p) == [0, 3, 4]
    @test get_bitmask(p, 1, UInt) == UInt(0b00111)
    @test get_bitmask(p, 2, UInt) == UInt(0b01000)

    @test typeof(compress(p, [1,1], UInt8)) === UInt8
    @test compress(p, [1,1]) == 0b1001
    @test compress(p, [2,1]) == 0b1010
    @test compress(p, [3,1]) == 0b1011
    # @test_throws ArgumentError compress(p, [0,0,0], UInt)
    # @test_throws ArgumentError compress(p, [-1,1])
    # @test_throws ArgumentError compress(p, [8,1])

    @testset "large" begin
        b1 = Boson(:b, 256)
        p1 = ParticleSector(b1, f)
        @test_throws ArgumentError compress(p1, [0,0], UInt8)
    end

    @test extract(p, 0b1010) == [2, 1]
    # @test_throws ArgumentError extract(p, 0b1111)
    @test extract(p, 0b10000000) == [0, 0]
end
