export ParticleState
export ParticleSite

export particlebitwidth

import ExactDiagonalization.bitwidth
import ExactDiagonalization.simplify

# Base.:(+)(x::F, y::F) where {F<:Fermion} = F(x.value ⊻ y.value)


struct ParticleState{
    P<:ParticleSector,      # P<:Tuple{Vararg{<:AbstractParticle}},
    QN<:Tuple{Vararg{<:AbstractQuantumNumber}},
}
    name::String
    particlesector::P
    quantumnumber::QN

    function ParticleState(
        name::AbstractString,
        particles::Tuple{Vararg{AbstractParticle}},
        quantumnumber::QN,
    ) where {QN}
        particlesector = ParticleSector(particles)
        P = typeof(particlesector)
        new{P, QN}(name, particlesector, quantumnumber)
    end
end

# particlebitwidth(::ParticleState{QN, P}) where {QN, P} = mapreduce(particlebitwidth, +, P.parameters)
# particlebitwidth(::Type{ParticleState{QN, P}}) where {QN, P} = mapreduce(particlebitwidth, +, P.parameters)

particlebitwidth(::ParticleState{P, QN}) where {P, QN} = particlebitwidth(P)
particlebitwidth(::Type{ParticleState{P, QN}}) where {P, QN} = particlebitwidth(P)

# particlecount(s::ParticleState{QN, <:NTuple{N, Any}}) where {QN, N} = N

# function particlebitwidthlist(state::ParticleState{QN, P}) where {QN, P}
#     return map(particlebitwidth, state.particles)
# end


struct ParticleSite{
    P<:ParticleSector,
    QN<:Tuple{Vararg{<:AbstractQuantumNumber}},
}
    states::Vector{ParticleState{P, QN}}
    occupationmap::Dict{P, Int}
    function ParticleSite(states::AbstractVector{ParticleState{P, QN}}) where {P, QN}
        occupationmap = Dict{P, Int}()
        for (i, s) in enumerate(states)
            if haskey(occupationmap, s.particlesector)
                throw(ArgumentError("duplicate particle count"))
            end
            occupationmap[s.particlesector] = i
        end
        return new{P, QN}(states, occupationmap)
    end
end

bitwidth(site::ParticleSite) = Int(ceil(log2(length(site.states))))
dimension(site::ParticleSite) = length(site.states)

export getstateindex
function getstateindex(site::ParticleSite{P, QN}, occupation::P) where {P<:ParticleSector, QN}
    return site.occupationmap[occupation]
end
