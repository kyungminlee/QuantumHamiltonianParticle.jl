using Test
using QuantumHamiltonianParticle

@testset "ladder particle" begin

    @testset "fermion" begin
        p = ParticleSector(Fermion(:f))
        site = ParticleSite([ParticleState(p, "_", [0]), ParticleState(p, "f", [1])])
        hs = ParticleHilbertSpace([site])
        hsr = represent(hs)
        dim = length(hsr.basis_list)

        for (c, cdag) in [
            (ParticleLadderUnit(p, 1, 1, ANNIHILATION),
             ParticleLadderUnit(p, 1, 1, CREATION)),
            (make_projector_operator(ParticleLadderUnit(p, 1, 1, ANNIHILATION)),
             make_projector_operator(ParticleLadderUnit(p, 1, 1, CREATION))),
        ]
            mat_c = zeros(Float64, (dim, dim))
            mat_cdag = zeros(Float64, (dim, dim))
            for (ivec, bvec) in enumerate(hsr.basis_list)
                for (bcol, ampl) in get_row_iterator(c, bvec)
                    icol = get(hsr.basis_lookup, bcol, -1)
                    (icol > 0) && (mat_c[ivec, icol] = ampl)
                end
                for (bcol, ampl) in get_row_iterator(cdag, bvec)
                    icol = get(hsr.basis_lookup, bcol, -1)
                    (icol > 0) && (mat_cdag[ivec, icol] = ampl)
                end
            end
            @test mat_c == [0.0 1.0; 0.0 0.0]
            @test mat_cdag == [0.0 0.0; 1.0 0.0]

            mat_c = zeros(Float64, (dim, dim))
            mat_cdag = zeros(Float64, (dim, dim))
            for (ivec, bvec) in enumerate(hsr.basis_list)
                for (brow, ampl) in get_column_iterator(c, bvec)
                    irow = get(hsr.basis_lookup, brow, -1)
                    (irow > 0) && (mat_c[irow, ivec] = ampl)
                end
                for (brow, ampl) in get_column_iterator(cdag, bvec)
                    irow = get(hsr.basis_lookup, brow, -1)
                    (irow > 0) && (mat_cdag[irow, ivec] = ampl)
                end
            end
            @test mat_c == [0.0 1.0; 0.0 0.0]
            @test mat_cdag == [0.0 0.0; 1.0 0.0]
        end
    end


    @testset "hardcore boson" begin
        p = ParticleSector(HardcoreBoson(:h))
        site = ParticleSite([ParticleState(p, "_", [0]), ParticleState(p, "f", [1])])
        hs = ParticleHilbertSpace([site])
        hsr = represent(hs)
        dim = length(hsr.basis_list)

        for (c, cdag) in [
            (ParticleLadderUnit(p, 1, 1, ANNIHILATION),
             ParticleLadderUnit(p, 1, 1, CREATION)),
            (make_projector_operator(ParticleLadderUnit(p, 1, 1, ANNIHILATION)),
             make_projector_operator(ParticleLadderUnit(p, 1, 1, CREATION))),
        ]
            mat_c = zeros(Float64, (dim, dim))
            mat_cdag = zeros(Float64, (dim, dim))
            for (ivec, bvec) in enumerate(hsr.basis_list)
                for (bcol, ampl) in get_row_iterator(c, bvec)
                    icol = get(hsr.basis_lookup, bcol, -1)
                    (icol > 0) && (mat_c[ivec, icol] = ampl)
                end
                for (bcol, ampl) in get_row_iterator(cdag, bvec)
                    icol = get(hsr.basis_lookup, bcol, -1)
                    (icol > 0) && (mat_cdag[ivec, icol] = ampl)
                end
            end
            @test mat_c == [0 1; 0 0]
            @test mat_cdag == [0 0; 1 0]

            mat_c = zeros(Float64, (dim, dim))
            mat_cdag = zeros(Float64, (dim, dim))
            for (ivec, bvec) in enumerate(hsr.basis_list)
                for (brow, ampl) in get_column_iterator(c, bvec)
                    irow = get(hsr.basis_lookup, brow, -1)
                    (irow > 0) && (mat_c[irow, ivec] = ampl)
                end
                for (brow, ampl) in get_column_iterator(cdag, bvec)
                    irow = get(hsr.basis_lookup, brow, -1)
                    (irow > 0) && (mat_cdag[irow, ivec] = ampl)
                end
            end
            @test mat_c == [0 1; 0 0]
            @test mat_cdag == [0 0; 1 0]
        end
    end


    @testset "boson" begin
        p = ParticleSector(Boson(:b, 3))
        site = ParticleSite([
            ParticleState(p, "0", [0]),
            ParticleState(p, "1", [1]),
            ParticleState(p, "2", [2]),
            ParticleState(p, "3", [3]),
        ])
        hs = ParticleHilbertSpace([site])
        hsr = represent(hs)
        dim = length(hsr.basis_list)

        for (c, cdag) in [
            (ParticleLadderUnit(p, 1, 1, ANNIHILATION),
             ParticleLadderUnit(p, 1, 1, CREATION)),
            (make_projector_operator(ParticleLadderUnit(p, 1, 1, ANNIHILATION)),
             make_projector_operator(ParticleLadderUnit(p, 1, 1, CREATION))),
        ]

            mat_c = zeros(Float64, (dim, dim))
            mat_cdag = zeros(Float64, (dim, dim))
            for (ivec, bvec) in enumerate(hsr.basis_list)
                for (bcol, ampl) in get_row_iterator(c, bvec)
                    icol = get(hsr.basis_lookup, bcol, -1)
                    (icol > 0) && (mat_c[ivec, icol] = ampl)
                end
                for (bcol, ampl) in get_row_iterator(cdag, bvec)
                    icol = get(hsr.basis_lookup, bcol, -1)
                    (icol > 0) && (mat_cdag[ivec, icol] = ampl)
                end
            end
            # https://en.wikipedia.org/wiki/Creation_and_annihilation_operators#Matrix_representation
            @test mat_c ≈ sqrt.([
                0 1 0 0;
                0 0 2 0;
                0 0 0 3;
                0 0 0 0
            ])
            @test mat_cdag ≈ sqrt.([
                0 0 0 0;
                1 0 0 0;
                0 2 0 0;
                0 0 3 0
            ])

            mat_c = zeros(Float64, (dim, dim))
            mat_cdag = zeros(Float64, (dim, dim))
            for (ivec, bvec) in enumerate(hsr.basis_list)
                for (brow, ampl) in get_column_iterator(c, bvec)
                    irow = get(hsr.basis_lookup, brow, -1)
                    (irow > 0) && (mat_c[irow, ivec] = ampl)
                end
                for (brow, ampl) in get_column_iterator(cdag, bvec)
                    irow = get(hsr.basis_lookup, brow, -1)
                    (irow > 0) && (mat_cdag[irow, ivec] = ampl)
                end
            end
            @test mat_c ≈ sqrt.([
                0 1 0 0;
                0 0 2 0;
                0 0 0 3;
                0 0 0 0
            ])
            @test mat_cdag ≈ sqrt.([
                0 0 0 0;
                1 0 0 0;
                0 2 0 0;
                0 0 3 0
            ])
        end
    end

    @testset "spin-1/2" begin
        p = ParticleSector(Spin(:s, 1))
        site = ParticleSite([
            ParticleState(p, "Sᶻ = 1/2", [0]),
            ParticleState(p, "Sᶻ =-1/2", [1]),
        ])
        hs = ParticleHilbertSpace([site])
        hsr = represent(hs)
        dim = length(hsr.basis_list)

        for (c, cdag) in [
            (ParticleLadderUnit(p, 1, 1, ANNIHILATION),
             ParticleLadderUnit(p, 1, 1, CREATION)),
            (make_projector_operator(ParticleLadderUnit(p, 1, 1, ANNIHILATION)),
             make_projector_operator(ParticleLadderUnit(p, 1, 1, CREATION))),
        ]
            mat_c = zeros(Float64, (dim, dim))
            mat_cdag = zeros(Float64, (dim, dim))
            for (ivec, bvec) in enumerate(hsr.basis_list)
                for (bcol, ampl) in get_row_iterator(c, bvec)
                    icol = get(hsr.basis_lookup, bcol, -1)
                    (icol > 0) && (mat_c[ivec, icol] = ampl)
                end
                for (bcol, ampl) in get_row_iterator(cdag, bvec)
                    icol = get(hsr.basis_lookup, bcol, -1)
                    (icol > 0) && (mat_cdag[ivec, icol] = ampl)
                end
            end
            # https://easyspin.org/easyspin/documentation/spinoperators.html
            @test mat_c ≈ sqrt.([
                0 1 ;
                0 0
            ])  # c is S+
            @test mat_cdag ≈ sqrt.([
                0 0;
                1 0
            ])  # cdag is S-

            mat_c = zeros(Float64, (dim, dim))
            mat_cdag = zeros(Float64, (dim, dim))
            for (ivec, bvec) in enumerate(hsr.basis_list)
                for (brow, ampl) in get_column_iterator(c, bvec)
                    irow = get(hsr.basis_lookup, brow, -1)
                    (irow > 0) && (mat_c[irow, ivec] = ampl)
                end
                for (brow, ampl) in get_column_iterator(cdag, bvec)
                    irow = get(hsr.basis_lookup, brow, -1)
                    (irow > 0) && (mat_cdag[irow, ivec] = ampl)
                end
            end
            @test mat_c ≈ sqrt.([
                0 1 ;
                0 0
            ])  # c is S+
            @test mat_cdag ≈ sqrt.([
                0 0;
                1 0
            ])  # cdag is S-

        end
    end


    @testset "spin-3/2" begin
        p = ParticleSector(Spin(:s, 3))
        site = ParticleSite([
            ParticleState(p, "Sᶻ = 3/2", [0]),
            ParticleState(p, "Sᶻ = 1/2", [1]),
            ParticleState(p, "Sᶻ =-1/2", [2]),
            ParticleState(p, "Sᶻ =-3/2", [3]),
        ])
        hs = ParticleHilbertSpace([site])
        hsr = represent(hs)
        dim = length(hsr.basis_list)

        for (c, cdag) in [
            (ParticleLadderUnit(p, 1, 1, ANNIHILATION),
             ParticleLadderUnit(p, 1, 1, CREATION)),
            (make_projector_operator(ParticleLadderUnit(p, 1, 1, ANNIHILATION)),
             make_projector_operator(ParticleLadderUnit(p, 1, 1, CREATION))),
        ]
            mat_c = zeros(Float64, (dim, dim))
            mat_cdag = zeros(Float64, (dim, dim))
            for (ivec, bvec) in enumerate(hsr.basis_list)
                for (bcol, ampl) in get_row_iterator(c, bvec)
                    icol = get(hsr.basis_lookup, bcol, -1)
                    (icol > 0) && (mat_c[ivec, icol] = ampl)
                end
                for (bcol, ampl) in get_row_iterator(cdag, bvec)
                    icol = get(hsr.basis_lookup, bcol, -1)
                    (icol > 0) && (mat_cdag[ivec, icol] = ampl)
                end
            end
            # https://easyspin.org/easyspin/documentation/spinoperators.html
            @test mat_c ≈ sqrt.([
                0 3 0 0;
                0 0 4 0;
                0 0 0 3;
                0 0 0 0
            ])  # c is S+
            @test mat_cdag ≈ sqrt.([
                0 0 0 0;
                3 0 0 0;
                0 4 0 0;
                0 0 3 0
            ])  # cdag is S-

            mat_c = zeros(Float64, (dim, dim))
            mat_cdag = zeros(Float64, (dim, dim))
            for (ivec, bvec) in enumerate(hsr.basis_list)
                for (brow, ampl) in get_column_iterator(c, bvec)
                    irow = get(hsr.basis_lookup, brow, -1)
                    (irow > 0) && (mat_c[irow, ivec] = ampl)
                end
                for (brow, ampl) in get_column_iterator(cdag, bvec)
                    irow = get(hsr.basis_lookup, brow, -1)
                    (irow > 0) && (mat_cdag[irow, ivec] = ampl)
                end
            end
            @test mat_c ≈ sqrt.([
                0 3 0 0;
                0 0 4 0;
                0 0 0 3;
                0 0 0 0
            ])  # c is S+
            @test mat_cdag ≈ sqrt.([
                0 0 0 0;
                3 0 0 0;
                0 4 0 0;
                0 0 3 0
            ])  # cdag is S-
        end
    end
end
