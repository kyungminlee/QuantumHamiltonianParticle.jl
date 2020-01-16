abstract type AbstractParticle end

struct Fermion{Species} <: AbstractParticle end
struct HardcoreBoson{Species} <: AbstractParticle end
struct Boson{Species, Max} <:AbstractParticle end

function pariclespecies(p::Type{F}) ::String where F<:AbstractParticle
  return p.parameters[1]
end

# Default is 1 for different types
exchangesign(::Type{T1}, ::Type{T2}) where {T1<:AbstractParticle, T2<:AbstractParticle} = 1
exchangesign(::Type{T}, ::Type{T}) where {T<:Fermion} = -1
exchangesign(::T1, ::T2) where {T1<:AbstractParticle, T2<:AbstractParticle} = exchangesign(T1, T2)

isfermion(::Type{T}) where {T<:AbstractParticle} = false
isfermion(::Type{T}) where {T<:Fermion} = true

maxoccupancy(::Type{<:Fermion}) = 1
maxoccupancy(::Type{<:HardcoreBoson}) = 1
maxoccupancy(::Type{<:Boson{Species, Max}}) where {Species, Max} = Max::Int
