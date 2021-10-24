export ParticleSector, ParticleIndex

import QuantumHamiltonian
import QuantumHamiltonian.bitwidth
import QuantumHamiltonian.compress
import QuantumHamiltonian.extract
import QuantumHamiltonian.get_bitmask

import QuantumHamiltonian.bitoffset
export numspecies, speciescount, getspecies, getspeciesname
export findspeciesindex

struct ParticleSector{
    P<:Tuple{Vararg{AbstractParticle}},
    # Cache Info
    E, # exchangesigns (tuple of integers)
    M, # maxoccupancies (tuple of integers)
    BW, # bitwidths (tuple of integers)
    BO, # bitoffsets (tuple of integers. cumsum of BW, starting with 0.)
}
    function ParticleSector(::Type{P}) where {P<:Tuple{Vararg{AbstractParticle}}}
        t = tuple(P.parameters...)
        M = maxoccupancy.(t)
        E = exchangesign.(t)
        BW = bitwidth.(t)
        BO = (0, cumsum(BW)...)
        return new{P, E, M, BW, BO}()
    end
    function ParticleSector(::P) where {P<:Tuple{Vararg{AbstractParticle}}}
        return ParticleSector(P)
    end
    function ParticleSector(p::Vararg{AbstractParticle})
        return ParticleSector(typeof(p))
    end
end


numspecies(::Type{<:ParticleSector{P}}) where {P} = tuplelength(P)
speciescount(::Type{<:ParticleSector{P}}) where {P} = tuplelength(P)
getspecies(::Type{P}) where {P<:ParticleSector} = tuple(P.parameters[1].parameters...)
getspecies(::Type{<:ParticleSector{P}}, index::Integer) where {P} = P.parameters[index]
getspeciesname(::Type{PS}, index::Integer) where {PS<:ParticleSector} = getspecies(PS, index).parameters[1]::Symbol

findspeciesindex(::P, args...) where {P<:ParticleSector} = findspeciesindex(P, args...)

# TODO: more efficient implementation
function findspeciesindex(::Type{<:ParticleSector{P}}, name::Symbol) where {P}
    i = findfirst(x -> getspeciesname(x) == name, P.parameters)
    return isnothing(i) ? -1 : Int(i)
end

# exchangesign(::Type{PS}, iptl::Integer) where {PS<:ParticleSector} = exchangesign(getspecies(PS, iptl))
exchangesign(::Type{<:ParticleSector{P, E}}, iptl::Integer) where {P, E} = E[iptl]
function exchangesign(::Type{PS}, iptl1::Integer, iptl2::Integer) where {PS<:ParticleSector}
    return iptl1 == iptl2 ? exchangesign(PS, iptl1) : 1
end

# occupation representaiton

# bitwidth(::Type{P}) where {P<:ParticleSector} = sum(bitwidth(p) for p in getspecies(P))::Int
# bitwidth(::Type{P}, iptl::Integer) where {P<:ParticleSector} = bitwidth(getspecies(P, iptl))::Int
bitwidth(::Type{PS}) where {PS<:ParticleSector} = sum(bitwidth(PS, iptl) for iptl in 1:numspecies(PS))::Int
bitwidth(::Type{<:ParticleSector{P, E, M, B}}, iptl::Integer) where {P, E, M, B} = B[iptl]::Int

# function bitoffset(::Type{P}, iptl::Integer)::Int where {P<:ParticleSector}
#     spec = getspecies(P)
#     offset = 0
#     for i in 1:(iptl-1)
#         offset += bitwidth(spec[i])
#     end
#     return offset
# end

function bitoffset(::Type{<:ParticleSector{P, E, M, BW, BO}}, iptl::Integer)::Int where {P, E, M, BW, BO}
    return BO[iptl]
end

function bitoffset(::Type{P}, iptl::Integer, isite::Integer)::Int where {P<:ParticleSector}
    bw = bitwidth(P)
    return bitoffset(P, iptl) + (bw * (isite-1))
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
    ::Type{BR}=UInt,
)::BR where {P<:ParticleSector, BR<:Unsigned}
  return make_bitmask(bitwidth(P), BR)
end


function get_bitmask(
    ::Type{P},
    iptl::Integer,
    ::Type{BR}=UInt,
)::BR where {P<:ParticleSector, BR<:Unsigned}
  offset = bitoffset(P, iptl)
  return make_bitmask(offset+bitwidth(P, iptl), offset, BR)
end


# function get_bitmask(
#     ::Type{P},
#     iptl::Integer,
#     isite::Integer,
#     ::Type{BR}=UInt,
# )::BR where {P<:ParticleSector, BR<:Unsigned}
#     offset = bitoffset(P, iptl)
#     bm = make_bitmask(offset+bitwidth(P, iptl), offset, BR)
#     bw = bitwidth(P)
#     return bm << (bw*(isite-1))
# end


"""
    get_bitmask(phs, [iptl, isite])

Get the bit mask for the particles `iptl` at sites `isite`.
`iptl` or `isite` can either be integer, a vector of integers, or colon `:`.
Bitwise or is taken over list of iptl.
"""
function get_bitmask(
    ::Type{PS},
    iptl::Integer,
    isite::Integer,
    ::Type{BR}=UInt,
)::BR where {PS<:ParticleSector, BR<:Unsigned}
    @boundscheck !(1<= iptl <= speciescount(PS)) && throw(BoundsError(PS.parameters, iptl))
    bm = get_bitmask(PS, iptl, BR)
    bw = bitwidth(PS)
    return bm << (bw * (isite-1))
end

function get_bitmask(
    ::Type{PS},
    iptls::AbstractVector{<:Integer},
    isites::AbstractVector{<:Integer},
    ::Type{BR}=UInt,
)::BR where {PS<:ParticleSector, BR<:Unsigned}
    @boundscheck for iptl in iptls
        !(1<= iptl <= speciescount(PS)) && throw(BoundsError(PS.parameters, iptl))
    end
    bm = mapreduce(iptl -> get_bitmask(PS, iptl, BR), |, iptls; init=zero(BR))
    bw = bitwidth(PS)
    return mapreduce(isite -> (bm << (bw * (isite-1))), |, isites; init=zero(BR))
end

function get_bitmask(
    ::Type{PS},
    iptl::Integer,
    isites::AbstractVector{<:Integer},
    ::Type{BR}=UInt,
)::BR where {PS<:ParticleSector, BR<:Unsigned}
    @boundscheck !(1<= iptl <= speciescount(PS)) && throw(BoundsError(PS.parameters, iptl))
    bm = get_bitmask(PS, iptl, BR)
    bw = bitwidth(PS)
    return mapreduce(isite -> (bm << (bw * (isite-1))), |, isites; init=zero(BR))
end

function get_bitmask(
    ::Type{PS},
    iptls::AbstractVector{<:Integer},
    isite::Integer,
    ::Type{BR}=UInt,
)::BR where {PS<:ParticleSector, BR<:Unsigned}
    @boundscheck for iptl in iptls
        !(1<= iptl <= speciescount(PS)) && throw(BoundsError(PS.parameters, iptl))
    end
    bm = mapreduce(iptl -> get_bitmask(PS, iptl, BR), |, iptls; init=zero(BR))
    bw = bitwidth(PS)
    return bm << (bw * (isite-1))
end


function get_parity_bitmask(
    ::Type{PS},
    iptl::Integer,
    isite::Integer,
    binary_type::Type{BR}=UInt,
) where {PS<:ParticleSector, BR<:Unsigned}
    if isfermion(getspecies(PS, iptl))
        bm_mask = zero(BR)
        bmp = get_bitmask(PS, iptl, BR)
        bw = bitwidth(PS)
        for jsite in 0:(isite-2)
            bm_mask |= bmp << (bw*jsite) 
        end
        return bm_mask
    else
        return zero(BR)
    end
end


function get_occupancy(
    ::Type{PS},
    iptl::Integer,
    isite::Integer,
    bvec::BR,
) where {PS<:ParticleSector, BR<:Unsigned}
    bm = get_bitmask(PS, iptl, BR)
    return Int( ((bvec >> (bitwidth(PS)*(isite-1))) & bm) >> bitoffset(PS, iptl) )
end


function get_occupancy(
    ::Type{PS},
    iptl::Integer,
    bvec::BR,
) where {PS<:ParticleSector, BR<:Unsigned}
    bm = get_bitmask(PS, iptl, BR)
    return Int( (bvec & bm) >> bitoffset(PS, iptl) )
end




"""
    set_occupancy(phs, iptl, isite, bvec::Unsigned, count)

Set occupancy of particle `iptl` at site `isite` for the given basis state `bvec` to `count`.
"""
function set_occupancy(
    ::Type{PS},
    iptl::Integer,
    isite::Integer,
    bvec::BR,
    count::Integer,
) where {PS, BR}
    @boundscheck !(0 <= count <= maxoccupancy(getspecies(PS, iptl))) && throw(ArgumentError("count out of bounds"))
    bw = bitwidth(PS)
    bm = get_bitmask(PS, iptl, BR) << (bw*(isite-1))
    bo = bitoffset(PS, iptl, isite)
    return (bvec & ~bm) | (BR(count) << bo)
end


function compress(
    ::Type{PS},
    occupancy::AbstractVector{<:Integer},
    ::Type{BR}=UInt,
)::BR where {PS<:ParticleSector, BR<:Unsigned}
    if length(occupancy) != speciescount(PS)
        throw(ArgumentError("length of occupancy vector should match the number of particles"))
    elseif sizeof(BR) * 8 < bitwidth(PS)
        throw(ArgumentError("type $BR is too short to represent the particle sector"))
    end

    out = zero(BR)
    offset = 0
    # for (i, (p, n)) in enumerate(zip(getspecies(PS), occupancy))
    for (p, n) in zip(getspecies(PS), occupancy)
        # n < 0 && throw(ArgumentError("occupancy should be non-negative"))
        # n > maxoccupancy(p) && throw(ArgumentError("occupancy ($n) should be no greater than the maxoccupancy of particle ($p)"))
        out |= BR(n) << offset
        offset += bitwidth(p)
    end
    return out
end


function extract(::Type{P}, occbin::BR)::Vector{Int} where {P<:ParticleSector, BR<:Unsigned}
    n_particles = speciescount(P)
    occ = Vector{Int}(undef, n_particles)
    for (i, p) in enumerate(getspecies(P))
        mask = (one(BR) << bitwidth(p)) - 1
        n = Int(occbin & mask)
        # @boundscheck if n > maxoccupancy(p)
        #     throw(ArgumentError("Occupation $n of $p exceeds $(maxoccupancy(p))"))
        # end
        occ[i] = n
        occbin >>= bitwidth(p)
    end
    return occ
end

function extract(::Type{P}, occbin::BR, occ::AbstractVector{<:Integer})::Vector{Int} where {P<:ParticleSector, BR<:Unsigned}
    for (i, p) in enumerate(getspecies(P))
        mask = (one(BR) << bitwidth(p)) - 1
        n = Int(occbin & mask)
        occ[i] = n
        occbin >>= bitwidth(p)
    end
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
    :bitwidth, :bitoffset, :get_bitmask, :get_parity_bitmask,
    :get_occupancy, :set_occupancy,
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
