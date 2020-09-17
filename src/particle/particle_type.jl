abstract type AbstractParticle end

export Fermion, HardcoreBoson, Boson, Spin
export species
export exchangesign
export isfermion, isboson
export maxoccupancy

import ExactDiagonalization.bitwidth
import ExactDiagonalization.compress
import ExactDiagonalization.extract
import ExactDiagonalization.get_bitmask

export bitoffset

struct Fermion{Species}<:AbstractParticle
    Fermion(species::Symbol) = new{species}()
    Fermion(species::AbstractString) = new{Symbol(species)}()
    Fermion{Species}() where {Species} = new{Species}()
end

struct HardcoreBoson{Species}<:AbstractParticle
    HardcoreBoson(species::Symbol) = new{species}()
    HardcoreBoson(species::AbstractString) = new{Symbol(species)}()
    HardcoreBoson{Species}() where {Species} = new{Species}()
end

struct Boson{Species, Max}<:AbstractParticle
    Boson(species::Symbol, max::Integer) = new{species, max}()
    Boson(species::AbstractString, max::Integer) = new{Symbol(species), max}()
    Boson{Species, Max}() where {Species, Max} = new{Species, Max}()
end

struct Spin{Species, N}<:AbstractParticle
    Spin(species::Symbol, max::Integer) = new{species, max}()
    Spin(species::AbstractString, max::Integer) = new{Symbol(species), max}()
    Spin{Species, Max}() where {Species, Max} = new{Species, Max}()
end


# getspecies(p::Type{<:AbstractParticle}) = p.parameters[1]

# in units of 1/q where q is the ``charge'' of the particle
# particleflux(::Type{<:Boson}) = 0//1
# particleflux(::Type{<:HardcoreBoson}) = 0//1
# particleflux(::Type{<:Fermion}) = 1//2 # TODO: is 1 a good number? or should i say 1/2 == pi/(2pi) ?

exchangesign(::Type{<:Boson}) = 1
exchangesign(::Type{<:HardcoreBoson}) = 1
exchangesign(::Type{<:Fermion}) = -1
exchangesign(::Type{<:Spin}) = 1

isfermion(::Type{<:AbstractParticle}) = false
isfermion(::Type{<:Fermion}) = true

isboson(::Type{<:AbstractParticle}) = false
isboson(::Type{<:HardcoreBoson}) = true
isboson(::Type{<:Boson}) = true
isboson(::Type{<:Spin}) = false

isspin(::Type{<:AbstractParticle}) = false
isspin(::Type{<:HardcoreBoson}) = false
isspin(::Type{<:Boson}) = false
isspin(::Type{<:Spin}) = true

maxoccupancy(::Type{<:Spin{S, M}}) where {S, M} = M
maxoccupancy(::Type{<:Boson{S, M}}) where {S, M} = M
maxoccupancy(::Type{<:HardcoreBoson}) = 1
maxoccupancy(::Type{<:Fermion}) = 1

# TODO: Edge case: when maxoccupancy is 0.
bitwidth(::Type{P}) where {P<:AbstractParticle} = Int(ceil(log2(maxoccupancy(P)+1)))

for fname in [:exchangesign, :isfermion, :isboson, :maxoccupancy, :bitwidth, :bitoffset]
    @eval begin
        ($fname)(::T, args...) where {T<:AbstractParticle} = ($fname)(T, args...)
    end
end
