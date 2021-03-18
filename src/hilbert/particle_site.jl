export ParticleState
export ParticleSite

using QuantumHamiltonian
import QuantumHamiltonian.bitwidth
import QuantumHamiltonian.qntype
import QuantumHamiltonian.dimension
import QuantumHamiltonian.bitwidth
import QuantumHamiltonian.get_state
import QuantumHamiltonian.get_state_index
import QuantumHamiltonian.get_quantum_number
import QuantumHamiltonian.quantum_number_sectors
import QuantumHamiltonian.compress


"""
    ParticleState{PS, BR, QN}

A state, represented by particle occupation.
"""
struct ParticleState{PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}
    name::String
    occupancy_binary::BR
    quantum_number::QN

    @doc """
        ParticleState(::Type{PS}, name, occvec, quantum_number, ::Type{BR}=UInt)

    Create a particle state

    # Arguments
    - `PS`: particle sector
    - `name`: name of the state
    - `occvec::AbstractVector{<:Integer}`: occupation Vector
    - `quantum_number`: quantum number (tuple or integer)
    - `BR`: binary type
    """
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

    @doc """
        ParticleState(::Type{PS}, name, occvec, ::Type{BR}=UInt)

    Create a particle state with no quantum number.

    # Arguments
    - `PS`: particle sector
    - `name`: name of the state
    - `occvec::AbstractVector{<:Integer}`: occupation Vector
    - `BR`: binary type
    """
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


# decorated with particle
"""
    ParticleSite{PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}

Particle site.

# Example
```
ParticleSite([state1, state2, state3, state4])
```
"""
struct ParticleSite{PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}
    states::Vector{ParticleState{PS, BR, QN}}
    state_lookup::Dict{BR, Int}

    function ParticleSite(states::AbstractVector{ParticleState{PS, BR, QN}}) where {QN, PS, BR}
        n_particles, n_states = numspecies(PS), length(states)
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

dimension(site::ParticleSite) = length(site.states)

function get_state_index(site::ParticleSite, binrep::Unsigned)
    return site.state_lookup[binrep]
end

function get_state(site::ParticleSite, binrep::Unsigned)
    return site.states[get_state_index(site, binrep)]
end

function quantum_number_sectors(site::ParticleSite{PS, BR, QN})::Vector{QN} where {PS, BR, QN}
    return sort(collect(Set([state.quantum_number for state in site.states])))
end

function get_quantum_number(site::ParticleSite{PS, BR, QN}, state_index::Integer)::QN where {PS, BR, QN}
    return site.states[state_index].quantum_number
end


for fname in [
    :exchangesign,
    :numspecies, :speciescount, :getspecies, :getspeciesname,
]
    @eval begin
        ($fname)(::Type{ParticleState{PS, BR, QN}}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
        ($fname)(::Type{ParticleSite{PS, BR, QN}}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
        ($fname)(::ParticleState{PS, BR, QN}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
        ($fname)(::ParticleSite{PS, BR, QN}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
    end
end

for fname in [
    :bitwidth, :bitoffset, :get_bitmask,
    :compress, :extract,
]
    @eval begin
        ($fname)(::Type{ParticleState{PS, BR, QN}}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
        ($fname)(::Type{ParticleSite{PS, BR, QN}}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
        ($fname)(::ParticleState{PS, BR, QN}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
        ($fname)(::ParticleSite{PS, BR, QN}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
    end
end
