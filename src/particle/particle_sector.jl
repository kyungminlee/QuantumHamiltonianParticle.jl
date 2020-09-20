export ParticleSector, ParticleIndex

import ExactDiagonalization.bitwidth
import ExactDiagonalization.compress
import ExactDiagonalization.extract
import ExactDiagonalization.get_bitmask

export bitoffset, bitoffset
export numspecies, speciescount, getspecies, getspeciesname


struct ParticleSector{P<:Tuple{Vararg{AbstractParticle}}}
    function ParticleSector(::Type{P}) where {P<:Tuple{Vararg{AbstractParticle}}}
        return new{P}()
    end
    function ParticleSector(::P) where {P<:Tuple{Vararg{AbstractParticle}}}
        return new{P}()
    end
    function ParticleSector(p::Vararg{AbstractParticle})
        return new{typeof(p)}()
    end
end


numspecies(::Type{ParticleSector{P}}) where {P} = tuplelength(P)
speciescount(::Type{ParticleSector{P}}) where {P} = tuplelength(P)
getspecies(::Type{P}) where {P<:ParticleSector} = tuple(P.parameters[1].parameters...)
getspecies(::Type{ParticleSector{P}}, index::Integer) where {P} = P.parameters[index]
getspeciesname(::Type{ParticleSector{P}}, index::Integer) where {P} = getspecies(ParticleSector{P}, index).parameters[1]::Symbol

exchangesign(::Type{PS}, iptl::Integer) where {PS<:ParticleSector} = exchangesign(getspecies(PS, iptl))
function exchangesign(::Type{PS}, iptl1::Integer, iptl2::Integer) where {PS<:ParticleSector}
    return iptl1 == iptl2 ? exchangesign(PS, iptl1) : 1
end


# occupation representaiton

bitwidth(::Type{P}) where {P<:ParticleSector} = sum(bitwidth(p) for p in getspecies(P))
bitwidth(::Type{P}, iptl::Integer) where {P<:ParticleSector} = bitwidth(getspecies(P, iptl))


function bitoffset(::Type{P}, iptl::Integer)::Int where {P<:ParticleSector}
    spec = getspecies(P)
    offset = 0
    for i in 1:(iptl-1)
        offset += bitwidth(spec[i])
    end
    return offset
end


function bitoffset(::Type{P}) where {P<:ParticleSector}
    spec = getspecies(P)
    nptls = length(spec)
    offset = 0
    out = Vector{Int}(undef, length(spec)+1)
    for i in 1:nptls
        out[i] = offset
        offset += bitwidth(spec[i])
    end
    out[end] = offset
    return out
end


function get_bitmask(
    ::Type{P},
    iptl::Integer,
    binary_type::Type{BR}=UInt,
)::BR where {P<:ParticleSector, BR<:Unsigned}
  offset = bitoffset(P, iptl)
  return make_bitmask(offset+bitwidth(P, iptl), offset, BR)
end


function compress(
    ::Type{P},
    occupancy::AbstractVector{<:Integer},
    binary_type::Type{BR}=UInt,
)::BR where {P<:ParticleSector, BR<:Unsigned}
    if length(occupancy) != speciescount(P)
        throw(ArgumentError("length of occupancy vector should match the number of particles"))
    elseif sizeof(BR) * 8 < bitwidth(P)
        throw(ArgumentError("type $BR is too short to represent the particle sector"))
    end

    out = zero(BR)
    offset = 0
    for (i, (p, n)) in enumerate(zip(getspecies(P), occupancy))
        n < 0 && throw(ArgumentError("occupancy should be non-negative"))
        n > maxoccupancy(p) && throw(ArgumentError("occupancy ($n) should be no greater than the maxoccupancy of particle ($p)"))
        out |= BR(n) << offset
        offset += bitwidth(p)
    end
    return out
end


function extract(::Type{P}, occbin::BR)::Vector{Int} where {P<:ParticleSector, BR<:Unsigned}
    offset = 0
    n_particles = speciescount(P)
    occ = Vector{Int}(undef, n_particles)
    for (i, p) in enumerate(getspecies(P))
        mask = (one(BR) << bitwidth(p)) - 1
        n = Int(occbin & mask)
        @boundscheck if n > maxoccupancy(p)
            throw(ArgumentError("Occupation $n of $p exceeds $(maxoccupancy(p))"))
        end
        occ[i] = n
        occbin >>= bitwidth(p)
    end
    # @assert(iszero(occbin))
    return occ
end



for fname in [
    :exchangesign,
    :numspecies, :speciescount, :getspecies, :getspeciesname,
]
    @eval begin
        ($fname)(p::P, args...) where {P<:ParticleSector} = ($fname)(P, args...)
    end
end

for fname in [
    :bitwidth, :bitoffset, :get_bitmask,
    :compress, :extract,
]
    @eval begin
        ($fname)(p::P, args...) where {P<:ParticleSector} = ($fname)(P, args...)
    end
end


# # Mainly for the ladder operators
# struct ParticleIndex{P<:ParticleSector}
#     index::Int

#     function ParticleIndex(::Type{P}, index::Integer) where {P<:ParticleSector}
#         NP = num_particle_species(P)
#         (0 < index <= NP) || throw(ArgumentError("index should be 0 < index <= NP"))
#         new{P}(index)
#     end

#     ParticleIndex(::P, args...) where {P<:ParticleSector} = ParticleIndex(P, args...)
# end


# Base.convert(::Type{T}, pi::ParticleIndex) where {T<:Integer} = convert(T, pi.index)


# Base.isless(lhs::ParticleIndex{P}, rhs::ParticleIndex{P}) where P = Base.isless(lhs.index, rhs.index)
# Base.:(==)(lhs::ParticleIndex{P}, rhs::ParticleIndex{P}) where P = lhs.index == rhs.index


# getspecies(p::ParticleIndex{P}) where {P<:ParticleSector} = P.parameters[1].parameters[p.index]
# getspeciesname(p::ParticleIndex{P}) where {P<:ParticleSector} = getspecies(p).parameters[1]
# bitwidth(p::ParticleIndex{P}) where {P<:ParticleSector} = bitwidth(P, p.index)
# bitoffset(p::ParticleIndex{P}) where {P<:ParticleSector} = bitoffset(P, p.index)
# get_bitmask(p::ParticleIndex{P}) where {P<:ParticleSector} = get_bitmask(P, p.index)
