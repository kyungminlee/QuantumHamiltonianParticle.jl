using Test
using Particle
using ExactDiagonalization

@testset "test projection" begin
    p1 = ParticleProjectorUnitOperator(0b0101, 0b0001, 0b0100, 0b0010, 2)
    p2 = ParticleProjectorUnitOperator(0b0101, 0b0001, 0b0100, 0b0010, 3.0)
    p3 = ParticleProjectorUnitOperator(UInt(0b0101), UInt(0b0001), UInt(0b0100), UInt(0b0010), 3.0)

    l1 = [p1, 3.0]
    @test eltype(l1) == ParticleProjectorUnitOperator{UInt8, Float64}
    @test l1[1] == p1
    @test l1[2] == ParticleProjectorUnitOperator(0b0000, 0b0000, 0b0000, 0b0000, 3.0)

    l2 = [p1, NullOperator()]
    @test eltype(l2) == ParticleProjectorUnitOperator{UInt8, Int}
    @test l2[1] == p1
    @test l2[2] == ParticleProjectorUnitOperator(0b0000, 0b0000, 0b0000, 0b0000, 0)

    l3 = [p1, p2]
    @test eltype(l3) == ParticleProjectorUnitOperator{UInt8, Float64}
    @test l3[1] == p1
    @test l3[2] == p2

    l4 = [p1, p3]
    @test eltype(l4) == ParticleProjectorUnitOperator{UInt, Float64}
    @test l4[1] == p1
    @test l4[2] == p3

    l5 = [p1, p1 + p1]
    @test eltype(l5) == ParticleProjectorSumOperator{UInt8, Int}
    @test l5[1] == ParticleProjectorSumOperator([p1])
    @test l5[2] == ParticleProjectorSumOperator([p1, p1])

    l6 = [p3, p1 + p1]
    @test eltype(l6) == ParticleProjectorSumOperator{UInt, Float64}
    @test l6[1] == ParticleProjectorSumOperator([p3])
    @test l6[2] == ParticleProjectorSumOperator([p1, p1])

    l7 = [p1 + p1,  p3 + p3 + p3]
    @test eltype(l7) == ParticleProjectorSumOperator{UInt, Float64}
    @test l7[1] == ParticleProjectorSumOperator([
        ParticleProjectorUnitOperator(UInt(0b0101), UInt(0b0001), UInt(0b0100), UInt(0b0010), 2.0),
        ParticleProjectorUnitOperator(UInt(0b0101), UInt(0b0001), UInt(0b0100), UInt(0b0010), 2.0),
    ])
end
