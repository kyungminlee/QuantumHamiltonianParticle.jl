using ExactDiagonalization
using Particle
using Test
using LinearAlgebra
using Formatting

@testset "spinless fermion" begin
    fermion = Fermion{:f}()
    particle_sector = make_particle_sector(fermion)
    site = ParticleSite([
        ParticleState(particle_sector, "_", [0], (0,)),
        ParticleState(particle_sector, "c", [1], (1,)),
    ])
    nsites = 8
    hs = ParticleHilbertSpace([site for i in 1:nsites])

    c_dag(i) = LadderUnitOperator(1, i, CREATION)
    c(i) = LadderUnitOperator(1, i, ANNIHILATION)

    @testset "checks" begin
        @test length(hs.sites) == nsites
        @test num_particle_species(hs) == 1
        for isite in 1:nsites
            @test get_bitmask(hs, 1, isite) == one(UInt) << (isite-1)
        end
    end

    @testset "application" begin
        @testset "ladder unit" begin
            #             87654321
            bvec = UInt(0b01010110)
            out = get_column_iterator(hs, c(8), bvec)
            @test isempty(out)

            out = get_column_iterator(hs, c_dag(8), bvec)
            @test first(out) == (UInt(0b11010110)=>1.0)

            out = get_column_iterator(hs, c_dag(6), bvec)
            (newbvec, ampl) = first(out)
            @test first(out) == (UInt(0b01110110) => -1.0)

            out = get_column_iterator(hs, c_dag(4), bvec)
            @test first(out) == (UInt(0b01011110) => 1.0)
        end

        @testset "ladder product" begin
            #             87654321
            bvec = UInt(0b01010110)

            out = get_column_iterator(hs, c_dag(8)*c(7), bvec)
            @test first(out) == (UInt(0b10010110) => 1.0)

            out = get_column_iterator(hs, c_dag(8)*c(5), bvec)
            @test first(out) == (UInt(0b11000110) => -1.0)

            out = get_column_iterator(hs, c_dag(8)*c(4), bvec)
            @test isempty(out)

            out = get_column_iterator(hs, c_dag(8)*c(3), bvec)
            @test first(out) == (UInt(0b11010010) => 1.0)
        end

        @testset "ladder sum" begin
            #             87654321
            bvec = UInt(0b01010110)
            out = get_column_iterator(hs, c_dag(8)*c(5) + c_dag(8)*c(3) * 0.5, bvec)
            @test Dict(out) == Dict(UInt(0b11000110) => -1.0, UInt(0b11010010) => 0.5)
        end


        @testset "hopping" begin
            t1 = 1.0
            t2 = 0.2
            hopping1 = sum(
                let j = mod(i, nsites) + 1
                    k = mod(i+1, nsites) + 1
                    t1 * (c_dag(i)*c(j) + c_dag(j)*c(i)) + t2 * (c_dag(i)*c(k) + c_dag(k)*c(i))
                end
                    for i in 1:nsites
            )
            hopping2 = make_projection_operator(hs, hopping1)

            for bcol in UInt(0):UInt(1<<nsites-1)
                a = Dict( get_column_iterator(hs, hopping1, bcol) )
                b = Dict( get_column_iterator(hopping2, bcol) )
                @test a == b
            end

            @testset "hopping matrix" begin
                hopping_matrix0 = t1 * [
                    0 1 0 0 0 0 0 1;
                    1 0 1 0 0 0 0 0;
                    0 1 0 1 0 0 0 0;
                    0 0 1 0 1 0 0 0;
                    0 0 0 1 0 1 0 0;
                    0 0 0 0 1 0 1 0;
                    0 0 0 0 0 1 0 1;
                    1 0 0 0 0 0 1 0;
                ] + t2 * [
                    0 0 1 0 0 0 1 0;
                    0 0 0 1 0 0 0 1;
                    1 0 0 0 1 0 0 0;
                    0 1 0 0 0 1 0 0;
                    0 0 1 0 0 0 1 0;
                    0 0 0 1 0 0 0 1;
                    1 0 0 0 1 0 0 0;
                    0 1 0 0 0 1 0 0;
                ]

                basis_list = [UInt(0x1) << (i-1) for i in 1:nsites]
                basis_lookup = Dict(b => i for (i, b) in enumerate(basis_list))
                let
                    hopping_matrix1 = zeros(Float64, (nsites, nsites))
                    for (icol, bcol) in enumerate(basis_list)
                        for (brow, ampl) in get_column_iterator(hs, hopping1, bcol)
                            irow = basis_lookup[brow]
                            hopping_matrix1[irow, icol] += ampl
                        end
                    end
                    @test hopping_matrix1 == hopping_matrix0
                end

                let
                    hopping_matrix2 = zeros(Float64, (nsites, nsites))
                    for (icol, bcol) in enumerate(basis_list)
                        for (brow, ampl) in get_column_iterator(hopping2, bcol)
                            irow = basis_lookup[brow]
                            hopping_matrix2[irow, icol] += ampl
                        end
                    end
                    @test hopping_matrix2 == hopping_matrix0
                end
            end

            @testset "eigenspectrum: three particles" begin
                eigenenergies0 = [
                    2*t1*(cos(2π*i/nsites) + cos(2π*j/nsites) + cos(2π*k/nsites)) + 
                    2*t2*(cos(4π*i/nsites) + cos(4π*j/nsites) + cos(4π*k/nsites))
                        for i in 0:(nsites-1) for j in 0:(i-1) for k in 0:(j-1)
                ]
                sort!(eigenenergies0)

                basis_list = [
                    (UInt(0x1) << (i-1)) | (UInt(0x1) << (j-1)) | (UInt(0x1) << (k-1))
                        for i in 1:nsites for j in 1:(i-1) for k in 1:(j-1)
                ]
                basis_lookup = Dict(b => i for (i, b) in enumerate(basis_list))
                nbasis = length(basis_list)

                let
                    hopping_matrix1 = zeros(Float64, (nbasis, nbasis))
                    for (icol, bcol) in enumerate(basis_list)
                        for (brow, ampl) in get_column_iterator(hs, hopping1, bcol)
                            irow = basis_lookup[brow]
                            hopping_matrix1[irow, icol] += ampl
                        end
                    end
                    @test ishermitian(hopping_matrix1)
                    eigenenergies1 = eigvals(Hermitian(hopping_matrix1))
                    @test isapprox(eigenenergies0, eigenenergies1)
                end
                let
                    hopping_matrix2 = zeros(Float64, (nbasis, nbasis))
                    for (icol, bcol) in enumerate(basis_list)
                        for (brow, ampl) in get_column_iterator(hopping2, bcol)
                            irow = basis_lookup[brow]
                            hopping_matrix2[irow, icol] += ampl
                        end
                    end
                    @test ishermitian(hopping_matrix2)
                    eigenenergies2 = eigvals(Hermitian(hopping_matrix2))
                    @test isapprox(eigenenergies0, eigenenergies2)
                end
            end
        end

        #=
        t = 1.0

        hopping = sum(let j = mod(i, nsites)+1
            t * ( c_dag(j)*c(i) + c_dag(i) * c(j))
        end
            for i in 1:nsites)
        bvec = UInt(0x8)

        for i in 1:10
            out = apply(hs, hopping, bvec)
            println("# $i")
            for (bvec, ampl) in out
                printfmtln("{} => {}", bitstring(bvec), ampl)
            end
        end

        state = SparseState{Float64, UInt}(bvec => 1.0)
        println("# 0")
        for (b, a) in state.components
            printfmtln("{:8s} => {}", bitstring(b)[end-7:end], a)
        end

        for i in 1:4
            state = apply(hs, hopping, state)
            println("# $i")
            for (b, a) in state.components
                printfmtln("{:8s} => {}", bitstring(b)[end-7:end], a)
            end
        end
        =#

    end

    #=
    for i in 3:5
        println("# c†($i)")
        prettyprintln(make_projection_operator(hs, c_dag(i)))
        println("# c($i)")
        prettyprintln(make_projection_operator(hs, c(i)))
    end
    a = make_projection_operator(hs, c_dag(5))
    b = make_projection_operator(hs, c(2))
    prettyprintln(b*a)
    prettyprintln(a*b)
    prettyprintln(make_projection_operator(hs, c_dag(5) * c(2)))
    =#

end
