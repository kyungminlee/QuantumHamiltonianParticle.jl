export AbstractParticle
export Fermion
export ParticleSector

export particlebitwidth
export numparticletype
export particletypecount


abstract type AbstractParticle end

struct Fermion{Species}<:AbstractParticle
    value::Int
    function Fermion{Species}(value::Integer) where {Species}
        !isa(Species, Symbol) && throw(ArgumentError("Species needs to be a Symbol"))
        value < 0 && throw(ArgumentError("particle number has to be nonnegative"))
        return new{Species}(value)
    end
end

particlebitwidth(::Fermion) = 1
particlebitwidth(::Type{<:Fermion}) = 1

"""
    ParticleSector

Keeps track of particle count
"""
struct ParticleSector{P<:Tuple{Vararg{AbstractParticle}}}
    particlecounts::P
    function ParticleSector(count::Vararg{AbstractParticle})
        P = typeof(count)
        return new{P}(count)
    end
    function ParticleSector(count::Tuple{Vararg{AbstractParticle}})
        P = typeof(count)
        return new{P}(count)
    end
end

particlebitwidth(::PS) where {PS<:ParticleSector} = particlebitwidth(PS)
particlebitwidth(::Type{ParticleSector{P}}) where P = mapreduce(particlebitwidth, +, P.parameters)

function particlebitwidth(::Type{ParticleSector{Tuple{P1}}}) where {P1}
    particlebitwidth(P1)
end

function particlebitwidth(::Type{ParticleSector{Tuple{P1, P2}}}) where {P1, P2}
    particlebitwidth(P1) + particlebitwidth(P2)
end

function particlebitwidth(::Type{ParticleSector{Tuple{P1, P2, P3}}}) where {P1, P2, P3}
    particlebitwidth(P1) + particlebitwidth(P2) + particlebitwidth(P3)
end

particletypecount(::P) where {P<:ParticleSector} = particletypecount(P)
numparticletype(::P) where {P<:ParticleSector} = numparticletype(P)

particletypecount(::Type{<:ParticleSector{<:NTuple{N, Any}}}) where {N} = N
numparticletype(::Type{<:ParticleSector{<:NTuple{N, Any}}}) where {N} = N
