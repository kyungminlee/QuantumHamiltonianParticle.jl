using ExactDiagonalization
using Particle
using Test
using LinearAlgebra
using Formatting

@testset "spinless fermion" begin
    fermion = Fermion{:f}()
    particle_sector = make_particle_sector(fermion)

    c_dag(i) = LadderUnitOperator(1, i, CREATION)
    c(i) = LadderUnitOperator(1, i, ANNIHILATION)


    site = ParticleSite([
        ParticleState(particle_sector, "_", [0], (0,)),
        ParticleState(particle_sector, "f", [1], (1,)),
    ])

    nsites = 8
    hs = ParticleHilbertSpace([site for i in 1:nsites])

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
            (newbvec, ampl) = first(out)
            @test newbvec == UInt(0b11010110)
            @test ampl == 1

            out = get_column_iterator(hs, c_dag(6), bvec)
            (newbvec, ampl) = first(out)
            @test newbvec == UInt(0b01110110)
            @test ampl == -1
        end

        @testset "ladder product" begin
            #             87654321
            bvec = UInt(0b01010110)

            out = get_column_iterator(hs, c_dag(8)*c(7), bvec)
            (newbvec, ampl) = first(out)
            @test newbvec == UInt(0b10010110)
            @test ampl == 1

            out = get_column_iterator(hs, c_dag(8)*c(5), bvec)
            (newbvec, ampl) = first(out)
            @test newbvec == UInt(0b11000110)
            @test ampl == -1

            out = get_column_iterator(hs, c_dag(8)*c(4), bvec)
            @test isempty(out)

            out = get_column_iterator(hs, c_dag(8)*c(3), bvec)
            (newbvec, ampl) = first(out)
            @test newbvec == UInt(0b11010010)
            @test ampl == 1
        end

        @testset "ladder sum" begin
            #             87654321
            bvec = UInt(0b01010110)

            out = get_column_iterator(hs, c_dag(8)*c(5) + c_dag(8)*c(3) * 0.5, bvec)
            @test Dict(out) == Dict(UInt(0b11000110) => -1.0, UInt(0b11010010) => 0.5)
        end


        @testset "hopping" begin
            hopping = sum(
                let j = mod(i, nsites) + 1
                    c_dag(i)*c(j) + c_dag(j)*c(i)
                end
                    for i in 1:nsites
            )
            hopping2 = make_projection_operator(hs, hopping)

            for bcol in UInt(0):UInt(1<<nsites-1)
                a = Dict( get_column_iterator(hs, hopping, bcol) )
                b = Dict( get_column_iterator(hopping2, bcol) )
                @test a == b
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
