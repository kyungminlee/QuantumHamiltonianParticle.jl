using Test
using ExactDiagonalization
using LatticeTools
using Particle
using Random

@testset "ParticleLadderEmbedding" begin
    @testset "Constructor & Comparison" begin
    end

    @testset "Methods" begin
    end

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

        hsr = represent(hs)
        rng = MersenneTwister(0)
        for op in [cdag(1,2), cdag(2,2), c(1,2), c(2,2), cdag(1,3)*c(1,1), cdag(2,3) * c(2,1), cdag(2,3)*c(2,1)*2 + c(1,1)*0.5]
            for bvec in rand(rng, hsr.basis_list, 8)
                @test collect(get_column_iterator(hs, op, bvec)) == collect(get_column_iterator(embed(hs, op), bvec))
                @test collect(get_row_iterator(hs, op, bvec)) == collect(get_row_iterator(embed(hs, op), bvec))
            end
            for brow in rand(rng, hsr.basis_list, 3), bcol in rand(rng, hsr.basis_list, 3)
                @test get_element(hs, op, brow, bcol) == get_element(embed(hs, op), brow, bcol)
            end
        end
    end

    @testset "Unary Operation" begin
    end

    @testset "Binary Operation" begin
    end

    @testset "apply!" begin
    end

    @testset "simplify" begin
    end
end
