using ExactDiagonalization
import ExactDiagonalization.bitwidth

export ParticleState, ParticleSite


struct ParticleState{PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}
  name::String
  occupancy_binary::BR
  quantum_number::QN

  function ParticleState(::Type{PS},
                         name::AbstractString,
                         occupancy::AbstractVector{<:Integer},
                         quantum_number::QN,
                         binary_type::Type{BR}=UInt) where {
                            PS<:ParticleSector,
                            BR<:Unsigned,
                            QN<:Tuple{Vararg{<:AbstractQuantumNumber}},
                          }
    occ_bin = compress(PS, occupancy, BR)
    new{PS, BR, QN}(name, occ_bin, quantum_number)
  end

  function ParticleState(::Type{PS},
                         name::AbstractString,
                         occupancy::AbstractVector{<:Integer},
                         binary_type::Type{BR}=UInt) where {PS<:ParticleSector, BR<:Unsigned}
    return ParticleState(PS, name, occupancy, (), BR)
  end

  function ParticleState(::Type{PS},
                         name::AbstractString,
                         occupancy::AbstractVector{<:Integer},
                         quantum_number::Integer,
                         binary_type::Type{BR}=UInt) where {PS<:ParticleSector, BR<:Unsigned}
    return ParticleState(PS, name, occupancy, (quantum_number,), BR)
  end

  ParticleState(::PS, args...) where {PS<:ParticleSector} = ParticleState(PS, args...)
end




export ParticleSite

# decorated with particle
struct ParticleSite{PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}
  states::Vector{ParticleState{PS, BR, QN}}
  state_lookup::Dict{BR, Int}

  function ParticleSite(states::AbstractVector{ParticleState{PS, BR, QN}}) where {QN, PS, BR}
    n_states = length(states)
    n_particles = num_particle_species(PS)

    lookup = Dict{BR, Int}()
    for (i, s) in enumerate(states)
      occbin = s.occupancy_binary
      haskey(lookup, occbin) && throw(ArgumentError("Duplicate occupancy $(s.occbin)"))
      lookup[occbin] = i
    end
    return new{PS, BR, QN}(states, lookup)
  end
end

bitwidth(site::ParticleSite{PS, BR, QN}) where {PS, BR, QN} = bitwidth(PS)
