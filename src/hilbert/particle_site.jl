using ExactDiagonalization
import ExactDiagonalization.bitwidth

export ParticleState, ParticleSite

import ExactDiagonalization.qntype
import ExactDiagonalization.dimension
import ExactDiagonalization.bitwidth
import ExactDiagonalization.get_state
import ExactDiagonalization.get_quantum_number
import ExactDiagonalization.quantum_number_sectors
import ExactDiagonalization.compress


struct ParticleState{PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}
    name::String
    occupancy_binary::BR
    # occupancy::Vector{Int}
    quantum_number::QN

    function ParticleState(
        ::Type{PS},
        name::AbstractString,
        occvec::AbstractVector{<:Integer},
        quantum_number::QN,
        ::Type{BR}=UInt,
    ) where {
        PS<:ParticleSector,
        BR<:Unsigned,
        QN<:Tuple{Vararg{<:AbstractQuantumNumber}},
    }
        occbin = compress(PS, occvec, BR)
        new{PS, BR, QN}(name, occbin, quantum_number)
    end

    function ParticleState(
        ::Type{PS},
        name::AbstractString,
        occvec::AbstractVector{<:Integer},
        ::Type{BR}=UInt,
    ) where {PS<:ParticleSector, BR<:Unsigned}
        return ParticleState(PS, name, occvec, (), BR)
    end

    function ParticleState(
        ::Type{PS},
        name::AbstractString,
        occvec::AbstractVector{<:Integer},
        quantum_number::Integer,
        ::Type{BR}=UInt,
    ) where {PS<:ParticleSector, BR<:Unsigned}
        return ParticleState(PS, name, occvec, (quantum_number,), BR)
    end

    ParticleState(::PS, args...) where {PS<:ParticleSector} = ParticleState(PS, args...)
end

function Base.:(==)(lhs::ParticleState{PS, BR, QN}, rhs::ParticleState{PS, BR, QN}) where {PS, BR, QN}
    return lhs.name == rhs.name && lhs.occupancy_binary == rhs.occupancy_binary && lhs.quantum_number == rhs.quantum_number
end


qntype(::Type{ParticleState{PS, BR, QN}}) where {PS, BR, QN} = QN


export ParticleSite

# decorated with particle
struct ParticleSite{PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}
    states::Vector{ParticleState{PS, BR, QN}}
    state_lookup::Dict{BR, Int}

    function ParticleSite(states::AbstractVector{ParticleState{PS, BR, QN}}) where {QN, PS, BR}
        n_states = length(states)
        n_particles = speciescount(PS)

        lookup = Dict{BR, Int}()
        for (i, s) in enumerate(states)
            occbin = s.occupancy_binary
            haskey(lookup, occbin) && throw(ArgumentError("Duplicate occupancy $(s.occbin)"))
            lookup[occbin] = i
        end
        return new{PS, BR, QN}(states, lookup)
    end
end

qntype(::Type{ParticleSite{PS, BR, QN}}) where {PS, BR, QN} = QN

# speciescount(::P) where {P<:ParticleSite} = speciescount(P)
# speciescount(::Type{ParticleSite{PS, BR, QN}}) where {PS, BR, QN} = speciescount(PS)
# bitwidth(site::ParticleSite{PS, BR, QN}) where {PS, BR, QN} = bitwidth(PS)

dimension(site::ParticleSite) = length(site.states)

function get_state(site::ParticleSite, binrep::Unsigned)
    return site.states[get_state_index(site, binrep)]
end

function get_state_index(site::ParticleSite, binrep::Unsigned)
    return site.state_lookup[binrep]
end

function compress(site::ParticleSite{PS, BR, QN}, state_index::Integer) where {PS, BR, QN}
    return site.states[state_index]
end

function quantum_number_sectors(site::ParticleSite{PS, BR, QN})::Vector{QN} where {PS, BR, QN}
    return sort(collect(Set([state.quantum_number for state in site.states])))
end

function get_quantum_number(site::ParticleSite{PS, BR, QN}, state_index::Integer)::QN where {PS, BR, QN}
    return site.states[state_index].quantum_number
end


for fname in [
    :exchangesign,
    :bitwidth, :bitoffset, :get_bitmask,
    :compress, :extract,
    :numspecies, :speciescount, :getspecies, :getspeciesname,
]
    @eval begin
        ($fname)(p::Type{ParticleState{PS, BR, QN}}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
        ($fname)(p::Type{ParticleSite{PS, BR, QN}}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
    end
end
