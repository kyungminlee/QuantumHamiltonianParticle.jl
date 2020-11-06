using ExactDiagonalization
using Particle
using Test
using LinearAlgebra
using Formatting

fermion = Fermion{:f}()
particle_sector = make_particle_sector(fermion)

c_dag(i) = ParticleLadderUnit(1, i, CREATION)
c(i) = ParticleLadderUnit(1, i, ANNIHILATION)


site = ParticleSite([
    ParticleState(particle_sector, "_", [0], (0,)),
    ParticleState(particle_sector, "f", [1], (1,)),
])

nsites = 8
hs = ParticleHilbertSpace([site for i in 1:nsites])

#=
a = make_projector_operator(hs, c_dag(5))
b = make_projector_operator(hs, c(2))
b*a |> prettyprintln
a*b |> prettyprintln
make_projector_operator(hs, c_dag(5) * c(2)) |> prettyprintln
make_projector_operator(hs, c_dag(5) * c(2) + c_dag(2) * c(5)) |> prettyprintln
=#

#get_column_iterator(hs, c(1), UInt(0x1))
#get_column_iterator(hs, c(2), UInt(0x1))
hop = c_dag(3) * c(1) + c_dag(1) * c(3)
hop2 = make_projector_operator(hs, hop)

#get_column_iterator(hs, hop, UInt(0x1))

# Must Test

hopping = sum(
    let j = mod(i, nsites) + 1
        c_dag(i)*c(j) + c_dag(j)*c(i)
    end
        for i in 1:nsites
)
hopping2 = make_projector_operator(hs, hopping)

#@show collect( get_column_iterator(hs, hop, UInt(0x1)) )
#@show collect( get_column_iterator(hop2, UInt(0x1)) )
@testset "hopping comparison" begin
    for bcol in UInt(0):UInt(1<<8-1)
    #let bcol = UInt(0x8)
        #@show bitstring(bcol)
        a = Dict( get_column_iterator(hs, hopping, bcol) )
        b = Dict( get_column_iterator(hopping2, bcol) )
        # @show a
        # @show b
        @test a == b
    end
    #@show bitstring(0x5)

    #@show make_projector_operator(hs, c_dag(8)*c(1))

end
#get_column_iterator(hs, c_dag(3) * c(1) + c_dag(4) * c(1), UInt(0x1))
