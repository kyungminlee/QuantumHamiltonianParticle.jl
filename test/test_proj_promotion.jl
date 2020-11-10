using Test
using Particle

@testset "test projection" begin
    p1 = ParticleProjectorUnitOperator(0b0101, 0b0001, 0b0100, 0b0010, 2)
    p2 = ParticleProjectorUnitOperator(0b0101, 0b0001, 0b0100, 0b0010, 3.0)
    p3 = ParticleProjectorUnitOperator(UInt(0b0101), UInt(0b0001), UInt(0b0100), UInt(0b0010), 3.0)

    l0 = [p1, 3.0]
    @test typeof(l0) == Vector{ParticleProjectorUnitOperator{UInt8, Float64}}
    @test l0[1] == p1
    @test l0[2] == ParticleProjectorUnitOperator(0b0000, 0b0000, 0b0000, 0b0000, 3.0)

    l1 = [p1, p2]
    @test typeof(l1) == Vector{ParticleProjectorUnitOperator{UInt8, Float64}}

    l2 = [p1, p3]
    @test typeof(l2) == Vector{ParticleProjectorUnitOperator{UInt, Float64}}

    p4 = p1 + p1
    l3 = [p1, p4]
    @test typeof(l3) == Vector{ParticleProjectorSumOperator{UInt8, Int}}

    p4 = p1 + p1
    l4 = [p3, p4]
    @test typeof(l4) == Vector{ParticleProjectorSumOperator{UInt, Float64}}

    p5 = p3 + p3
    l5 = [p4, p5]
    @test typeof(l5) == Vector{ParticleProjectorSumOperator{UInt, Float64}}
end
