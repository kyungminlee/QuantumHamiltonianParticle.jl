using Test
using LinearAlgebra
using ExactDiagonalization
using Particle
using Random

@testset "Ladder Iterator" begin
    p = ParticleSector(Boson(:m, 2), Fermion(:f))
    cdag(i, j) = ParticleLadderUnit(p, i, j, CREATION)
    c(i, j)    = ParticleLadderUnit(p, i, j, ANNIHILATION)
    site = ParticleSite([
        ParticleState(p, "__", [0, 0], (0, 0)),
        ParticleState(p, "b_", [1, 0], (1, 0)),
        ParticleState(p, "B_", [2, 0], (2, 0)),
        ParticleState(p, "_f", [0, 1], (0, 1)),
        ParticleState(p, "bf", [1, 1], (1, 1)),
        ParticleState(p, "Bf", [2, 1], (2, 1)),
    ])
    hs = ParticleHilbertSpace([site, site, site])

    # ParticleLadderNull
    nop = ParticleLadderNull(p)
    @test isempty(collect(get_column_iterator(hs, nop, UInt(0b000_000_000))))
    @test isempty(collect(get_row_iterator(hs, nop, UInt(0b000_000_000))))

    # ParticleLadderUnit
    @test collect(get_column_iterator(hs, cdag(1,2), UInt(0b000_000_000))) == [UInt(0b000_001_000) => 1.0]
    @test collect(get_column_iterator(hs, cdag(1,2), UInt(0b000_000_001))) == [UInt(0b000_001_001) => 1.0]
    @test collect(get_column_iterator(hs, cdag(2,2), UInt(0b000_000_000))) == [UInt(0b000_100_000) => 1.0]
    @test collect(get_column_iterator(hs, cdag(2,2), UInt(0b000_000_001))) == [UInt(0b000_100_001) => 1.0]
    @test collect(get_column_iterator(hs, cdag(2,2), UInt(0b000_000_101))) == [UInt(0b000_100_101) => -1.0]

    @test collect(get_column_iterator(hs, c(1,2), UInt(0b000_001_000))) == [UInt(0b000_000_000) => 1.0]
    @test collect(get_column_iterator(hs, c(1,2), UInt(0b000_001_001))) == [UInt(0b000_000_001) => 1.0]
    @test collect(get_column_iterator(hs, c(2,2), UInt(0b000_100_000))) == [UInt(0b000_000_000) => 1.0]
    @test collect(get_column_iterator(hs, c(2,2), UInt(0b000_100_001))) == [UInt(0b000_000_001) => 1.0]
    @test collect(get_column_iterator(hs, c(2,2), UInt(0b000_100_101))) == [UInt(0b000_000_101) => -1.0]

    @test collect(get_row_iterator(hs, cdag(1,2), UInt(0b000_001_000))) == [UInt(0b000_000_000) => 1.0]
    @test collect(get_row_iterator(hs, cdag(1,2), UInt(0b000_001_001))) == [UInt(0b000_000_001) => 1.0]
    @test collect(get_row_iterator(hs, cdag(2,2), UInt(0b000_100_000))) == [UInt(0b000_000_000) => 1.0]
    @test collect(get_row_iterator(hs, cdag(2,2), UInt(0b000_100_001))) == [UInt(0b000_000_001) => 1.0]
    @test collect(get_row_iterator(hs, cdag(2,2), UInt(0b000_100_101))) == [UInt(0b000_000_101) => -1.0]

    @test collect(get_row_iterator(hs, c(1,2), UInt(0b000_000_000))) == [UInt(0b000_001_000) => 1.0]
    @test collect(get_row_iterator(hs, c(1,2), UInt(0b000_000_001))) == [UInt(0b000_001_001) => 1.0]
    @test collect(get_row_iterator(hs, c(2,2), UInt(0b000_000_000))) == [UInt(0b000_100_000) => 1.0]
    @test collect(get_row_iterator(hs, c(2,2), UInt(0b000_000_001))) == [UInt(0b000_100_001) => 1.0]
    @test collect(get_row_iterator(hs, c(2,2), UInt(0b000_000_101))) == [UInt(0b000_100_101) => -1.0]

    # ParticleLadderProduct
    @test isempty(collect(get_column_iterator(hs, cdag(1,3) * c(1,1), UInt(0b000_000_000))))
    @test collect(get_column_iterator(hs, cdag(1,3) * c(1,1), UInt(0b000_000_001))) == [UInt(0b001_000_000) => 1.0]
    @test collect(get_column_iterator(hs, cdag(1,3) * c(1,1), UInt(0b000_001_001))) == [UInt(0b001_001_000) => 1.0]

    @test isempty(collect(get_column_iterator(hs, cdag(2,3) * c(2,1), UInt(0b000_000_000))))
    @test collect(get_column_iterator(hs, cdag(2,3) * c(2,1), UInt(0b000_000_100))) == [UInt(0b100_000_000) => 1.0]
    @test collect(get_column_iterator(hs, cdag(2,3) * c(2,1), UInt(0b000_001_100))) == [UInt(0b100_001_000) => 1.0]
    @test collect(get_column_iterator(hs, cdag(2,3) * c(2,1), UInt(0b000_101_100))) == [UInt(0b100_101_000) => -1.0]
end

@testset "Iterator Comparison" begin
    p = ParticleSector(Boson(:m, 2), Fermion(:f), Spin("S=1", 2))
    cdag(i, j) = ParticleLadderUnit(p, i, j, CREATION)
    c(i, j)    = ParticleLadderUnit(p, i, j, ANNIHILATION)
    site = ParticleSite([
        ParticleState(p, "__↑", [0, 0, 0], (0, 0, 1)),
        ParticleState(p, "b_↑", [1, 0, 0], (1, 0, 1)),
        ParticleState(p, "B_↑", [2, 0, 0], (2, 0, 1)),
        ParticleState(p, "_f↑", [0, 1, 0], (0, 1, 1)),
        ParticleState(p, "bf↑", [1, 1, 0], (1, 1, 1)),
        ParticleState(p, "Bf↑", [2, 1, 0], (2, 1, 1)),
        ParticleState(p, "__⋅", [0, 0, 1], (0, 0, 0)),
        ParticleState(p, "b_⋅", [1, 0, 1], (1, 0, 0)),
        ParticleState(p, "B_⋅", [2, 0, 1], (2, 0, 0)),
        ParticleState(p, "_f⋅", [0, 1, 1], (0, 1, 0)),
        ParticleState(p, "bf⋅", [1, 1, 1], (1, 1, 0)),
        ParticleState(p, "Bf⋅", [2, 1, 1], (2, 1, 0)),
        ParticleState(p, "__↓", [0, 0, 2], (0, 0,-1)),
        ParticleState(p, "b_↓", [1, 0, 2], (1, 0,-1)),
        ParticleState(p, "B_↓", [2, 0, 2], (2, 0,-1)),
        ParticleState(p, "_f↓", [0, 1, 2], (0, 1,-1)),
        ParticleState(p, "bf↓", [1, 1, 2], (1, 1,-1)),
        ParticleState(p, "Bf↓", [2, 1, 2], (2, 1,-1)),
    ])

    hs = ParticleHilbertSpace([site, site,])
    hsr = represent(hs)

    @testset "ParticleLatticeUnit" begin
        rng = MersenneTwister(0)
        for iptl1 in 1:3, cop1 in [c, cdag]
            op1 = cop1(iptl1, 2)
            mr = Dict{Tuple{UInt, UInt}, Float64}()
            mc = Dict{Tuple{UInt, UInt}, Float64}()
            me = Dict{Tuple{UInt, UInt}, Float64}()
            for brow in hsr.basis_list
                for (bcol, ampl) in get_row_iterator(hs, op1, brow)
                    mr[(brow, bcol)] = get(mr, (brow, bcol), zero(Float64)) + ampl
                end
            end
            for bcol in hsr.basis_list
                for (brow, ampl) in get_column_iterator(hs, op1, bcol)
                    mc[(brow, bcol)] = get(mc, (brow, bcol), zero(Float64)) + ampl
                end
            end
            for brow in hsr.basis_list, bcol in hsr.basis_list
                ampl = get_element(hs, op1, brow, bcol)
                if !iszero(ampl)
                    me[(brow, bcol)] = ampl
                end
            end
            choptol!(mr, 1E-8)
            choptol!(mc, 1E-8)
            choptol!(me, 1E-8)
            @test mr == mc
            @test me == mc
        end
    end

    @testset "ParticleLatticeProduct" begin
        rng = MersenneTwister(0)
        for iptl1 in 1:3, iptl2 in 1:3, cop1 in [c, cdag], cop2 in [c, cdag]
            op1 = cop1(iptl1, 1) * cop2(iptl2, 2)
            mr = Dict{Tuple{UInt, UInt}, Float64}()
            mc = Dict{Tuple{UInt, UInt}, Float64}()
            me = Dict{Tuple{UInt, UInt}, Float64}()
            for brow in hsr.basis_list
                for (bcol, ampl) in get_row_iterator(hs, op1, brow)
                    mr[(brow, bcol)] = get(mr, (brow, bcol), zero(Float64)) + ampl
                end
            end
            for bcol in hsr.basis_list
                for (brow, ampl) in get_column_iterator(hs, op1, bcol)
                    mc[(brow, bcol)] = get(mc, (brow, bcol), zero(Float64)) + ampl
                end
            end
            for brow in hsr.basis_list, bcol in hsr.basis_list
                ampl = get_element(hs, op1, brow, bcol)
                if !iszero(ampl)
                    me[(brow, bcol)] = ampl
                end
            end
            choptol!(mr, 1E-8)
            choptol!(mc, 1E-8)
            choptol!(me, 1E-8)
            @test mr == mc
            @test me == mc
        end
    end

    @testset "ParticleLatticeSum" begin
        rng = MersenneTwister(0)
        for iptl1 in 1:3, iptl2 in 1:3, cop1 in [c, cdag], cop2 in [c, cdag]
            op1 = cop1(iptl1, 1) * cop2(iptl2, 2) * 0.5 - c(1, 2) * 2
            mr = Dict{Tuple{UInt, UInt}, Float64}()
            mc = Dict{Tuple{UInt, UInt}, Float64}()
            me = Dict{Tuple{UInt, UInt}, Float64}()
            for brow in hsr.basis_list
                for (bcol, ampl) in get_row_iterator(hs, op1, brow)
                    mr[(brow, bcol)] = get(mr, (brow, bcol), zero(Float64)) + ampl
                end
            end
            for bcol in hsr.basis_list
                for (brow, ampl) in get_column_iterator(hs, op1, bcol)
                    mc[(brow, bcol)] = get(mc, (brow, bcol), zero(Float64)) + ampl
                end
            end
            for brow in hsr.basis_list, bcol in hsr.basis_list
                ampl = get_element(hs, op1, brow, bcol)
                if !iszero(ampl)
                    me[(brow, bcol)] = ampl
                end
            end
            choptol!(mr, 1E-8)
            choptol!(mc, 1E-8)
            choptol!(me, 1E-8)
            @test mr == mc
            @test me == mc
        end
    end

    function test_product_iterators()
        rng = MersenneTwister(0)
        for bvec in rand(rng, hsr.basis_list, 5)
            for iptl1 in 1:3, iptl2 in 1:3, cop1 in [c, cdag], cop2 in [c, cdag]
                op1 = cop1(iptl1, 1)
                op2 = cop2(iptl2, 2)
                out1 = Dict{UInt, Float64}(collect(get_column_iterator(hs, op1 * op2, bvec))...)
                out2 = Dict{UInt, Float64}()
                for (bcol, ampl) in get_column_iterator(hs, op2, bvec)
                    for (brow, ampl2) in get_column_iterator(hs, op1, bcol)
                        out2[brow] = get(out2, brow, zero(Float64)) + ampl * ampl2
                    end
                end
                choptol!(out1, 1E-8)
                choptol!(out2, 1E-8)
                out1 != out2 && return false
            end
        end
        return true
    end
    @test test_product_iterators()

    function test_sum_iterators()
        rng = MersenneTwister(0)
        for bvec in rand(rng, hsr.basis_list, 5)
            for iptl1 in 1:3, iptl2 in 1:3, iptl3 in 1:3, cop1 in [c, cdag], cop2 in [c, cdag], cop3 in [c, cdag]
                op1 = cop1(iptl1, 1) * cop2(iptl2, 2)
                op2 = cop3(iptl2, 2)
                out1 = Dict{UInt, Float64}(collect(get_column_iterator(hs, op1 + op2, bvec))...)
                out2 = Dict{UInt, Float64}(collect(get_column_iterator(hs, op1, bvec))...)
                for (brow, ampl) in get_column_iterator(hs, op2, bvec)
                    out2[brow] = get(out2, brow, 0) + ampl
                end
                choptol!(out1, 1E-8)
                choptol!(out2, 1E-8)
                out1 != out2 && return false
            end
        end
        return true
    end
    @test test_sum_iterators()
end
