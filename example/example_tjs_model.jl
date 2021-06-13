using QuantumHamiltonian
using QuantumHamiltonianParticle
using Test
using LinearAlgebra


electronUp = Fermion{Symbol("e↑")}()
electronDn = Fermion{Symbol("e↓")}()
magnon = Boson{Symbol("m↓"), 2}()  # spin 1

particle_sector = ParticleSector(magnon, electronUp, electronDn)

# magnon_index = ParticleIndex(particle_sector, 1)
# electron_up_index = ParticleIndex(particle_sector, 2)
# electron_dn_index = ParticleIndex(particle_sector, 3)

# c_up_dag(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(electron_up_index, orbital, CREATION))
# c_up(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(electron_up_index, orbital, ANNIHILATION))
# c_dn_dag(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(electron_dn_index, orbital, CREATION))
# c_dn(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(electron_dn_index, orbital, ANNIHILATION))

# c_up_dag(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(2, orbital, CREATION))
# c_up(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(2, orbital, ANNIHILATION))
# c_dn_dag(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(3, orbital, CREATION))
# c_dn(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(3, orbital, ANNIHILATION))

c_up_dag(orbital::OrbitalType) where OrbitalType = ParticleLadderUnit(particle_sector, 2, orbital, CREATION)
c_up(orbital::OrbitalType) where OrbitalType = ParticleLadderUnit(particle_sector, 2, orbital, ANNIHILATION)
c_dn_dag(orbital::OrbitalType) where OrbitalType = ParticleLadderUnit(particle_sector, 3, orbital, CREATION)
c_dn(orbital::OrbitalType) where OrbitalType = ParticleLadderUnit(particle_sector, 3, orbital, ANNIHILATION)

# sf_site = ParticleSite([
#     ParticleState(particle_sector, "↑,e↑", [0,1,0], 2+1),
#     ParticleState(particle_sector, "↑,e↓", [0,0,1], 2-1),
#     ParticleState(particle_sector, "↓,e↑", [1,1,0], -2+1),
#     ParticleState(particle_sector, "↓,e↓", [1,0,1], -2-1),
# ])

tjs_site = ParticleSite([   # charge and spin
    ParticleState(particle_sector, "↑,em", [0, 0, 0], (0, 1 + 0)),
    ParticleState(particle_sector, "↑,e↑", [0, 1, 0], (1, 1 + 1)),
    ParticleState(particle_sector, "↑,e↓", [0, 0, 1], (1, 1 - 1)),

    ParticleState(particle_sector, "0,em", [1, 0, 0], (0, 0 + 0)),
    ParticleState(particle_sector, "0,e↑", [1, 1, 0], (1, 0 + 1)),
    ParticleState(particle_sector, "0,e↓", [1, 0, 1], (1, 0 - 1)),

    ParticleState(particle_sector, "↓,em", [2, 0, 0], (0, -1 + 0)),
    ParticleState(particle_sector, "↓,e↑", [2, 1, 0], (1, -1 + 1)),
    ParticleState(particle_sector, "↓,e↓", [2, 0, 1], (1, -1 - 1)),
])

spin_site = ParticleSite([
    ParticleState(particle_sector, "↑", [0, 0, 0], (0, 1)),
    ParticleState(particle_sector, "↓", [1, 0, 0], (0, -1)),
])

@show tjs_site

hs = ParticleHilbertSpace([tjs_site, spin_site, tjs_site, tjs_site, tjs_site])

nptls = speciescount(hs)
nsites = length(hs.sites)

for isite in 1:nsites, iptl in 1:nptls
    print("$iptl $isite ")
    print(bitstring(get_bitmask(hs, iptl, isite)))
    println()
end


# function getfermionparity(hs::ParticleHilbertSpace, op::ParticleLadderUnit{<:Integer, <:Integer}, bvec::Unsigned)
#     bm_species = get_bitmask(hs, op.particle_index, :)
#     bm_site = get_bitmask(hs, op.particle_index, op.orbital)
#     bm_mask = ( bm_site - 1 ) & bm_species  # σᶻ in jordan wigner string
#     bm_parity = bm_mask & bvec
#     return count_ones(bm_parity) % 2
# end


bvec = UInt( 0b0000_1_1_00_1_1_00_0_0_00_1_1_10 )
#                   <--4-> <--3-> <--2-> <--1->

println("# Annihilation")
(newbvec, ampl) = apply(hs, c_up(3), bvec)
println(bvec |> bitstring)
println(newbvec |> bitstring)
println(ampl)

println("# Creation")
@show c_up_dag(4)
(newbvec, ampl) = apply(hs, c_up_dag(4), bvec)
println(bvec |> bitstring)
println(newbvec |> bitstring)
println(ampl)

println("# Density")
n3u = ParticleLadderProduct([c_up_dag(3), c_up(3)])
newbvec, ampl = apply(hs, n3u, bvec)
println(bvec |> bitstring)
println(newbvec |> bitstring)
println(ampl)


println("# Hopping 4->5")
hop = ParticleLadderProduct([c_up_dag(5), c_up(4)])
newbvec, ampl = apply(hs, hop, bvec)
println(bvec |> bitstring)
println(newbvec |> bitstring)
println(ampl)


println("# Hopping 4->2")
hop = ParticleLadderProduct([c_up_dag(2), c_up(4)])
newbvec, ampl = apply(hs, hop, bvec)
println(bvec |> bitstring)
println(newbvec |> bitstring)
println(ampl)
