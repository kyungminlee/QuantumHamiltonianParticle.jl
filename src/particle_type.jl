abstract type AbstractParticle end

export Fermion, HardcoreBoson, Boson
export particle_species
export exchangesign
export isfermion
export maxoccupancy

import ExactDiagonalization.bitwidth
import ExactDiagonalization.compress
import ExactDiagonalization.extract
import ExactDiagonalization.get_bitmask

export bitoffset

struct Fermion{Species}<:AbstractParticle end
struct HardcoreBoson{Species}<:AbstractParticle end
struct Boson{Species, Max}<:AbstractParticle end

particle_species(p::Type{<:AbstractParticle}) = p.parameters[1]

exchangesign(p::Type{<:Boson}) = 1
exchangesign(p::Type{<:HardcoreBoson}) = 1
exchangesign(p::Type{<:Fermion}) = -1

isfermion(p::Type{<:AbstractParticle}) = false
isfermion(p::Type{<:Fermion}) = true

isboson(::Type{<:AbstractParticle}) = false
isboson(::Type{<:HardcoreBoson}) = true
isboson(::Type{<:Boson}) = true

maxoccupancy(::Type{<:Boson{S, M}}) where {S, M} = M
maxoccupancy(::Type{<:HardcoreBoson}) = 1
maxoccupancy(::Type{<:Fermion}) = 1

# TODO: Edge case: when maxoccupancy is 0.
bitwidth(::Type{P}) where {P<:AbstractParticle} = Int(ceil(log2(maxoccupancy(P)+1)))

for fname in [:particle_species, :exchangesign, :isfermion, :isboson, :maxoccupancy, :bitwidth]
  @eval begin
    ($fname)(::T) where {T<:AbstractParticle} = ($fname)(T)
  end
end
