using ExactDiagonalization
using Particle
using Test
using LinearAlgebra

@testset "Particle" begin
    electronUp = Fermion{Symbol("e↑")}()
    electronDn = Fermion{Symbol("e↓")}()
    magnon = Boson{Symbol("m↓"), 2}()

    particle_sector = make_particle_sector(magnon, electronUp, electronDn)

    @show particle_sector
    @show typeof(particle_sector)
    @show sizeof(particle_sector)
    @show isbits(particle_sector)

    magnon_index = ParticleIndex(particle_sector, 1)
    electron_up_index = ParticleIndex(particle_sector, 2)
    electron_dn_index = ParticleIndex(particle_sector, 3)

    @show magnon_index
    @show electron_up_index
    @show isless(magnon_index, electron_up_index)
    @show isless(electron_up_index, magnon_index)

    c_up_dag(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(electron_up_index, orbital, CREATION))
    c_up(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(electron_up_index, orbital, ANNIHILATION))
    c_dn_dag(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(electron_dn_index, orbital, CREATION))
    c_dn(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(ParticleLadderUnit(electron_dn_index, orbital, ANNIHILATION))

    @testset "operator" begin

        @show c_up_dag("A")
        @show c_up("A")

        n_up = c_up_dag("A") * c_up("A")
        @show n_up
        prettyprint(n_up)
        println()
    end


    sf_site = ParticleSite([ParticleState(particle_sector, "↑,e↑", [0,1,0], 2+1),
                            ParticleState(particle_sector, "↑,e↓", [0,0,1], 2-1),
                            ParticleState(particle_sector, "↓,e↑", [2,1,0], -2+1),
                            ParticleState(particle_sector, "↓,e↓", [2,0,1], -2-1)])
    n_sites = 5
    phs = ParticleHilbertSpace([sf_site, sf_site,sf_site,sf_site,sf_site])
    pad = bitwidth(phs)
    @testset "hilbert" begin
        @show compress(particle_sector, [0,1,1])
        for i in 0:8
        @show i, extract(particle_sector, UInt(i))
        end
        for i in 1:length(sf_site.states)
        @show string(state2occbin(sf_site, i), base=2, pad=bitwidth(sf_site))
        @show state2occvec(sf_site, i)
        end

        @show bitwidth(phs)
        for i_site in 1:n_sites
        @show string(get_bitmask(phs, :, i_site), base=2, pad=bitwidth(phs))
        end

        @show string(get_bitmask(phs, 1:2, 1:3), base=2, pad=bitwidth(phs))

        @show string(statevec2occbin(phs, [1,1,1,1,1]), base=2, pad=bitwidth(phs))
        @show string(statevec2occbin(phs, [1,1,1,2,1]), base=2, pad=bitwidth(phs))
        @show string(statevec2occbin(phs, [1,1,1,3,1]), base=2, pad=bitwidth(phs))
        @show string(statevec2occbin(phs, [1,1,1,4,1]), base=2, pad=bitwidth(phs))

        @show statevec2occmat(phs, [1,1,1,1,1])
        @show statevec2occmat(phs, [1,1,1,2,1])
        @show statevec2occmat(phs, [1,1,1,3,1])
        @show statevec2occmat(phs, [1,1,1,4,1])


        @show occbin2statevec(phs, 0b0100_0100_0100_0100_0100)
        @show occbin2statevec(phs, 0b0100_1000_0100_0100_0100)
        @show occbin2statevec(phs, 0b0100_0110_0100_0100_0100)
        @show occbin2statevec(phs, 0b0100_1010_0100_0100_0100)

        @show statevec2locvec(phs, [1,1,1,1,1])
        @show statevec2locvec(phs, [1,1,1,2,1])
        @show statevec2locvec(phs, [1,1,1,3,1])
        @show statevec2locvec(phs, [1,1,1,4,1])

        @show locvec2statevec(phs, statevec2locvec(phs, [1,1,1,1,1]))
        @show locvec2statevec(phs, statevec2locvec(phs, [1,1,1,2,1]))
        @show locvec2statevec(phs, statevec2locvec(phs, [1,1,1,3,1]))
        @show locvec2statevec(phs, statevec2locvec(phs, [1,1,1,4,1]))

        @show get_bitmask(particle_sector, 1)
        @show get_bitmask(particle_sector, 2)
        @show get_bitmask(particle_sector, 3)

        @show string(get_bitmask(phs, 1, :), base=2, pad=pad)
        @show string(get_bitmask(phs, 2, :), base=2, pad=pad)
        @show string(get_bitmask(phs, 3, :), base=2, pad=pad)

        @show string(get_bitmask(phs, :, :), base=2, pad=pad)
    end

    @testset "moreop" begin
        cdag3 = represent(phs, ParticleLadderUnit(electron_up_index, 3, CREATION))
        c3 = represent(phs, ParticleLadderUnit(electron_up_index, 3, ANNIHILATION))

        println("cdag3")
        prettyprintln(cdag3)

        println("c3")
        prettyprintln(c3)

        n3 = cdag3 * c3

        println("n3")
        prettyprintln(n3)

        cup3dag = represent(phs, ParticleLadderUnit(electron_up_index, 3, CREATION))
        cdn3 = represent(phs, ParticleLadderUnit(electron_dn_index, 3, ANNIHILATION))

        println("cup3dag")
        prettyprintln(cup3dag)

        println("cdn3")
        prettyprintln(cdn3)
        flip = cup3dag * cdn3

        println("flip")
        prettyprintln(flip)

        println("make")
        prettyprintln(represent(phs, c_up_dag(3) * c_dn(3)))

        prettyprint(c_dn(3))
        println()
        prettyprint(adjoint(c_dn(3)))
        println()

        flip = c_up_dag(3) * c_dn(3)
        prettyprint(flip)
        println()
        prettyprint(adjoint(flip))
        println()

        #@show ishermitian(flip)
        flip2 = c_up_dag(3) * c_dn(3) + c_dn_dag(3) * c_up(3)
        prettyprint(flip2)
        println()
        prettyprint(adjoint(flip2))
        println()

        prettyprint((flip2 - adjoint(flip2)))
        println()
        prettyprint(simplify(flip2 - adjoint(flip2)))
        println()

        @show ishermitian(flip)
        @show ishermitian(flip2)

        flip2_rep = represent(phs, flip2)
        prettyprintln(flip2_rep)
    end
end




# #@show phs
#
# @show particle_decompose(phs, [1,2,3,4,1])
# foo = particle_decompose(phs, [1,2,3,4,1])
# @show particle_compose(phs, foo)
#
# @show particle_compose(phs, ([4,3], [1,3,5], [2,4]))
# @show particle_compose(phs, ([4,3], [1,5,3], [2,4]))
# @show particle_compose(phs, ([4,3], [5,3,1], [2,4]))
#
# @show particle_occupancy(phs, [1,1,1,1,1])
# @show particle_occupancy(phs, [1,2,3,3,4])
#
#






#=

P = Electron
cre(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(CreationOperator(P, orbital))
ann(orbital::OrbitalType) where OrbitalType = ParticleLadderSum(AnnihilationOperator(P, orbital))

hopping = sum(cre(mod(i, 4)+1) * ann(i) + cre(i) * ann(mod(i,4)+1) for i in 1:4)
#hopping = sum(cre(i) * cre(mod(i, 4)+1) + ann(mod(i,4)+1) * ann(i) for i in 1:4)

c_hopping = simplify(hopping)

println("hopping")
hopping |> prettyprint
println()
println()

println("c_hopping")
c_hopping |> prettyprint
println()
println()


h2 = cre(2) * ann(1) + cre(1) * ann(2)
for i in 2:4
    global h2
    h2 = h2 + cre(mod(i, 4)+1) * ann(i) + cre(i) * ann(mod(i,4)+1)
end

sample = cre(4) + cre(3) * 2.0 + cre(2) * (3im) + cre(1) * (4.0 * im)

#   Tuple{Vararg{Tuple{<:ParticleLadderProduct, <:Number}}}

new_sample = unify_type(sample)



cre(4) * cre(3) * cre(2) |> prettyprint
println()

cre(4) * cre(3) * cre(2) |> normal_order |> prettyprint
println()


normal_order(ann(3) * cre(3)) |> prettyprint
println()

(ann(3) * cre(3)) |> prettyprint
println()


ElectronUp = Fermion{:ElectronUp}
ElectronDn = Fermion{:ElectronDn}

tj_site = Site([State("Em"), State("Up"), State("Dn")])
tj_hilbert = HilbertSpace([tj_site, tj_site, tj_site])

spm = StateParticleMap(Tuple{ElectronUp, ElectronDn}, [(0,0), (1,0), (0,1)])

n_sites = 3


tj_particle_hilbert = ParticleHilbertSpace(tj_hilbert, [spm, spm, spm])

for indexarray in Iterators.product([1:3 for i in 1:n_sites]...)
  indexarray = collect(indexarray)
  p = particle_decompose(tj_particle_hilbert, indexarray)
  @show indexarray, p
end


=#










# AbstractQuantumNumber = Integer

# struct State{QN<:AbstractQuantumNumber}
#   name ::String
#   quantum_number ::QN

#   State(name ::AbstractString) = new{Int}(name, 0)
#   State(name ::AbstractString, quantum_number ::QN) where {QN} = new{QN}(name, quantum_number)
#   State{QN}(name ::AbstractString, quantum_number ::QN) where {QN} = new{QN}(name, quantum_number)
# end


# abstract type AbstractSite{QN} end

# struct ParticleSite{P<:AbstractParticle, QN} <: AbstractSite{QN}
#   states ::Vector{State{QN}}
#   function ParticleSite(::Type{P}, states::AbstractVector{State{QN}}) where {P<:AbstractParticle, QN}
#     if length(states) != maxoccupancy(P)+1
#       throw(ArgumentError("number of states $(length(states)) does not match the maximum occupancy of the particle $P $(maxoccupancy(P))"))
#     end
#     new{P, QN}(states)
#   end
# end

# getparticle_species(p::ParticleSite{P}) where P = P
# dimension(p::ParticleSite{P}) where P = maxoccupancy(P)+1
# bitwidth(p::ParticleSite) = Int(ceil(log2(dimension(p))))
# get_state(p::ParticleSite, binrep::U) where {U<:Unsigned} = Int(binrep)+1

# e = ParticleSite(Electron, [State("0e", SVector(0,0,0)), State("1e", SVector(1,0,0))])
# n = ParticleSite(Fermion{:Neutron}, [State("0n", SVector(0,0,0)), State("1n", SVector(0,1,0))])
# b = ParticleSite(Alpha, [State("0a", SVector(0,0,0)), State("1a", SVector(0,0,1))])

# struct HilbertSpace{QN} #<: AbstractHilbertSpace
#   sites ::Vector{AbstractSite{QN}}
#   bitwidths ::Vector{Int}
#   bitoffsets ::Vector{Int}

#   function HilbertSpace(sites::Vararg{<:AbstractSite{QN}}) where QN
#     sites = AbstractSite{QN}[s for s in sites]
#     bitwidths = map(bitwidth, sites)
#     bitoffsets = Int[0, cumsum(bitwidths)...]
#     # particle_types = []
#     # for s in sites
#     #   p = getparticle_species(s)
#     #   if ! (p in particle_types)
#     #     push!(particle_types, p)
#     #   end
#     # end
#     new{QN}(sites, bitwidths, bitoffsets)
#   end
# end

# h = HilbertSpace(e,e,n,e,e,n,b)

# # --

# #using ExactDiagonalization



# struct ParticleBasis{O<:Tuple}
#   basis ::O
#   ParticleBasis(o::O) where O = new{O}(o)
# end

# function prettyprint(bvec::ParticleBasis)
#   print("|")
#   print(join( basis, ","))
#   print("⟩")
# end






#-


;
