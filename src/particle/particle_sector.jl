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


@inline _Particles(::Type{<:ParticleSector{P, E, M, BW, BO}}) where {P, E, M, BW, BO} = P
@inline _ExchangeSigns(::Type{<:ParticleSector{P, E, M, BW, BO}}) where {P, E, M, BW, BO} = E
@inline _MaxOccupancies(::Type{<:ParticleSector{P, E, M, BW, BO}}) where {P, E, M, BW, BO} = M
@inline _BitWidths(::Type{<:ParticleSector{P, E, M, BW, BO}}) where {P, E, M, BW, BO} = BW
@inline _BitOffsets(::Type{<:ParticleSector{P, E, M, BW, BO}}) where {P, E, M, BW, BO} = BO


numspecies(::Type{PS}) where {PS<:ParticleSector} = tuplelength(_Particles(PS))
speciescount(::Type{PS}) where {PS<:ParticleSector} = tuplelength(_Particles(PS))
getspecies(::Type{PS}) where {PS<:ParticleSector} = tuple(_Particles(PS).parameters...)
getspecies(::Type{PS}, index::Integer) where {PS<:ParticleSector} = _Particles(PS).parameters[index]
getspeciesname(::Type{PS}, index::Integer) where {PS<:ParticleSector} = getspecies(PS, index).parameters[1]::Symbol

findspeciesindex(::PS, args...) where {PS<:ParticleSector} = findspeciesindex(PS, args...)

function findspeciesindex(::Type{PS}, name::Symbol) where {PS<:ParticleSector}
    # TODO: more efficient implementation
    i = findfirst(x -> getspeciesname(x) == name, _Particles(PS).parameters)
    return isnothing(i) ? -1 : Int(i)
end


exchangesign(::Type{PS}, iptl::Integer) where {PS<:ParticleSector} = _ExchangeSigns(PS)[iptl]::Int

function exchangesign(::Type{PS}, iptl1::Integer, iptl2::Integer) where {PS<:ParticleSector}
    return iptl1 == iptl2 ? exchangesign(PS, iptl1) : 1
end

maxoccupancy(::Type{PS}, iptl::Integer) where {PS<:ParticleSector} = _MaxOccupancies(PS)[iptl]::Int

# occupation representation

bitwidth(::Type{PS}) where {PS<:ParticleSector} = sum(bitwidth(PS, iptl) for iptl in 1:numspecies(PS))::Int

bitwidth(::Type{PS}, iptl::Integer) where {PS<:ParticleSector} = _BitWidths(PS)[iptl]::Int


bitoffset(::Type{PS}, iptl::Integer) where {PS<:ParticleSector} = _BitOffsets(PS)[iptl]::Int

function bitoffset(::Type{PS}, iptl::Integer, isite::Integer)::Int where {PS<:ParticleSector}
    return bitoffset(PS, iptl) + (bitwidth(PS) * (isite-1))
end

function bitoffset(::Type{PS}) where {PS<:ParticleSector}
    species = getspecies(PS)
    nptls = length(species)
    offset = 0
    out = Vector{Int}(undef, length(species)+1)
    for i in 1:nptls
        out[i] = offset
        offset += bitwidth(species[i])
    end
    out[end] = offset
    return out
end


function get_bitmask(::Type{PS}, ::Type{BR}=UInt)::BR where {PS<:ParticleSector, BR<:Unsigned}
    return make_bitmask(bitwidth(PS), BR)
end

function get_bitmask(::Type{PS}, iptl::Integer, ::Type{BR}=UInt)::BR where {PS<:ParticleSector, BR<:Unsigned}
    offset = bitoffset(PS, iptl)
    return make_bitmask(offset + bitwidth(PS, iptl), offset, BR)
end

# function get_bitmask(
#     ::Type{PS},
#     iptl::Integer,
#     isite::Integer,
#     ::Type{BR}=UInt,
# )::BR where {PS<:ParticleSector, BR<:Unsigned}
#     offset = bitoffset(PS, iptl)
#     bm = make_bitmask(offset + bitwidth(PS, iptl), offset, BR)
#     sbw = bitwidth(PS)
#     return bm << (sbw*(isite-1))
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
    sbw = bitwidth(PS)
    return bm << (sbw * (isite-1))
end

function get_bitmask(
    ::Type{PS},
    iptls::AbstractVector{<:Integer},
    isites::AbstractVector{<:Integer},
    ::Type{BR}=UInt,
)::BR where {PS<:ParticleSector, BR<:Unsigned}
    nptls = speciescount(PS)
    @boundscheck for iptl in iptls
        !(1 <= iptl <= nptls) && throw(BoundsError(PS.parameters, iptl))
    end
    @inbounds bm = mapreduce(iptl -> get_bitmask(PS, iptl, BR), |, iptls; init=zero(BR))
    sbw = bitwidth(PS)
    return mapreduce(isite -> (bm << (sbw * (isite-1))), |, isites; init=zero(BR))
end

function get_bitmask(
    ::Type{PS},
    iptl::Integer,
    isites::AbstractVector{<:Integer},
    ::Type{BR}=UInt,
)::BR where {PS<:ParticleSector, BR<:Unsigned}
    @boundscheck !(1<= iptl <= speciescount(PS)) && throw(BoundsError(PS.parameters, iptl))
    bm = get_bitmask(PS, iptl, BR)
    sbw = bitwidth(PS)
    return mapreduce(isite -> (bm << (sbw * (isite-1))), |, isites; init=zero(BR))
end

function get_bitmask(
    ::Type{PS},
    iptls::AbstractVector{<:Integer},
    isite::Integer,
    ::Type{BR}=UInt,
)::BR where {PS<:ParticleSector, BR<:Unsigned}
    nptls = speciescount(PS)
    @boundscheck for iptl in iptls
        !(1 <= iptl <= nptls) && throw(BoundsError(PS.parameters, iptl))
    end
    bm = mapreduce(iptl -> get_bitmask(PS, iptl, BR), |, iptls; init=zero(BR))
    sbw = bitwidth(PS)
    return bm << (sbw * (isite-1))
end


function get_parity_bitmask(
    ::Type{PS},
    iptl::Integer,
    isite::Integer,
    ::Type{BR}=UInt,
) where {PS<:ParticleSector, BR<:Unsigned}
    if isfermion(PS, iptl)
        bm_mask = zero(BR)
        bm = get_bitmask(PS, iptl, BR)
        sbw = bitwidth(PS)
        for jsite in 0:(isite-2)
            bm_mask |= bm << (sbw*jsite) 
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
    sbw = bitwidth(PS)
    bo = bitoffset(PS, iptl)
    return Int(((bvec >> (sbw*(isite-1))) & bm) >> bo)
end

function get_occupancy(
    ::Type{PS},
    iptl::Integer,
    bvec::BR,
) where {PS<:ParticleSector, BR<:Unsigned}
    bm = get_bitmask(PS, iptl, BR)
    bo = bitoffset(PS, iptl)
    return Int((bvec & bm) >> bo)
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
    @boundscheck !(0 <= count <= maxoccupancy(PS, iptl)) && throw(ArgumentError("count out of bounds"))
    sbw = bitwidth(PS)
    bm = get_bitmask(PS, iptl, BR) << (sbw*(isite-1))
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
    for (p, n) in zip(getspecies(PS), occupancy)
        # n < 0 && throw(ArgumentError("occupancy should be non-negative"))
        # n > maxoccupancy(p) && throw(ArgumentError("occupancy ($n) should be no greater than the maxoccupancy of particle ($p)"))
        out |= BR(n) << offset
        offset += bitwidth(p)
    end
    return out
end


function extract(::Type{PS}, occbin::BR)::Vector{Int} where {PS<:ParticleSector, BR<:Unsigned}
    nptls = speciescount(PS)
    occ = Vector{Int}(undef, nptls)    
    species = getspecies(PS)
    for (i, p) in enumerate(species)
        mask = (one(BR) << bitwidth(p)) - 1
        n = Int(occbin & mask)
        # @boundscheck if n > maxoccupancy(p) && throw(ArgumentError("Occupation $n of $p exceeds $(maxoccupancy(p))"))
        @inbounds occ[i] = n
        occbin >>= bitwidth(p)
    end
    return occ
end


function extract(::Type{PS}, occbin::BR, occ::AbstractVector{<:Integer}) where {PS<:ParticleSector, BR<:Unsigned}
    @warn "extract(PS, occbin, occ) deprecated. Use extract!(occ, PS, occbin) instead." maxlog=1
    species = getspecies(PS)
    @boundscheck length(occ) == length(species) || throw(ArgumentError("length(occ) $(length(occ)) != number of species $(length(species))"))
    for (i, p) in enumerate(species)
        mask = (one(BR) << bitwidth(p)) - 1
        n = Int(occbin & mask)
        @inbounds occ[i] = n
        occbin >>= bitwidth(p)
    end
    return occ
end


function extract!(occ::AbstractVector{<:Integer}, ::Type{PS}, occbin::BR) where {PS<:ParticleSector, BR<:Unsigned}
    species = getspecies(PS)
    @boundscheck length(occ) == length(species) || throw(ArgumentError("length(occ) $(length(occ)) != number of species $(length(species))"))
    for (i, p) in enumerate(species)
        mask = (one(BR) << bitwidth(p)) - 1
        n = Int(occbin & mask)
        @inbounds occ[i] = n
        occbin >>= bitwidth(p)
    end
    return occ
end


for fname in [
    :exchangesign,
    :numspecies, :speciescount, :getspecies, :getspeciesname,
]
    @eval begin
        ($fname)(::PS, args...) where {PS<:ParticleSector} = ($fname)(PS, args...)
    end
end

for fname in [
    :bitwidth, :bitoffset, :get_bitmask, :get_parity_bitmask,
    :get_occupancy, :set_occupancy,
    :compress, :extract,
]
    @eval begin
        ($fname)(::PS, args...) where {PS<:ParticleSector} = ($fname)(PS, args...)
    end
end

for fname in [
    :isfermion, :isboson, :isspin
]
    @eval begin
        ($fname)(::Type{PS}, i::Integer) where {PS<:ParticleSector} = ($fname)(getspecies(PS, i))
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
