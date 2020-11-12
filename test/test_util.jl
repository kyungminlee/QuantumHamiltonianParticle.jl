using Test
using Particle

@testset "util" begin
    @testset "tuple" begin
        P = Particle
        t1 = (1.0, 2, 3.0 + 4im)
        T1 = typeof(t1)
        @test P.tuplelength(T1)  === 3
        @test P.tupleone(T1)  === (1.0, 1, 1.0 + 0.0im)
        @test P.tupleone(T1)  ==  (1.0, 1, 1.0 + 0.0im)
        @test P.tuplezero(T1) === (0.0, 0, 0.0 + 0.0im)
        @test P.tuplezero(T1) ==  (0.0, 0, 0.0 + 0.0im)
        @test P.tupleone(T1)  !== (1, 1, 1)
        @test P.tupleone(T1)  ==  (1, 1, 1)
        @test P.tuplezero(T1) !== (0, 0, 0)
        @test P.tuplezero(T1) ==  (0, 0, 0)

        @test P.tuplelength(t1)  === 3
        @test P.tupleone(t1)  === (1.0, 1, 1.0 + 0.0im)
        @test P.tupleone(t1)  ==  (1.0, 1, 1.0 + 0.0im)
        @test P.tuplezero(t1) === (0.0, 0, 0.0 + 0.0im)
        @test P.tuplezero(t1) ==  (0.0, 0, 0.0 + 0.0im)
        @test P.tupleone(t1)  !== (1, 1, 1)
        @test P.tupleone(t1)  ==  (1, 1, 1)
        @test P.tuplezero(t1) !== (0, 0, 0)
        @test P.tuplezero(t1) ==  (0, 0, 0)

        @test P.tupleadd((1.0, 2, 3.0 + 4im), (5.0, 6, 7.0 + 8.0im)) === (6.0, 8, 10.0 + 12.0im)
        @test P.tuplesubtract((1.0, 2, 3.0 + 4im), (5.0, 6, 7.0 + 8.0im)) === (-4.0, -4, -4.0 - 4.0im)
    end
end
