using ExactDiagonalization
import ExactDiagonalization.get_bitmask
export get_bitmask_particle

export ParticleHilbertSpace

export get_fermion_parity
export get_occupancy, set_occupancy

import ExactDiagonalization.scalartype
# import ExactDiagonalization.valtype
import ExactDiagonalization.qntype
import ExactDiagonalization.basespace
import ExactDiagonalization.bitwidth
import ExactDiagonalization.bitoffset
import ExactDiagonalization.get_bitmask
import ExactDiagonalization.quantum_number_sectors
import ExactDiagonalization.get_quantum_number
#import ExactDiagonalization.

import ExactDiagonalization.AbstractHilbertSpace


# Add a decoration to the existing Hilbert space
struct ParticleHilbertSpace{PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}<:AbstractHilbertSpace
    sites::Vector{ParticleSite{PS, BR, QN}}
    bitwidths::Vector{Int}
    bitoffsets::Vector{Int}

    function ParticleHilbertSpace(sites::AbstractVector{ParticleSite{PS, BR, QN}}) where {PS, BR, QN}
        bitwidths = map(bitwidth, sites)
        bitoffsets = Int[0, cumsum(bitwidths)...]
        if sizeof(BR) * 8 < bitoffsets[end]
            throw(ArgumentError("type $BR too small to represent the hilbert space (need $(bitoffsets[end]) bits)"))
        end
        return new{PS, BR, QN}(sites, bitwidths, bitoffsets)
    end
end


function Base.:(==)(lhs::ParticleHilbertSpace, rhs::ParticleHilbertSpace)
    return lhs.sites == rhs.sites
end



scalartype(::Type{<:ParticleHilbertSpace}) = Bool
scalartype(::ParticleHilbertSpace) = Bool


Base.valtype(::Type{<:ParticleHilbertSpace}) = Bool
Base.valtype(::ParticleHilbertSpace) = Bool


qntype(::Type{ParticleHilbertSpace{PS, BR, QN}}) where {PS, BR, QN} = QN
qntype(::ParticleHilbertSpace{PS, BR, QN}) where {PS, BR, QN} = QN


basespace(hs::ParticleHilbertSpace) = hs


bitwidth(hs::ParticleHilbertSpace) = hs.bitoffsets[end]

# numspecies(::P) where {P<:ParticleHilbertSpace} = numspecies(P)
# numspecies(::Type{ParticleHilbertSpace{PS, BR, QN}}) where {PS, BR, QN} = numspecies(PS)
# speciescount(::P) where {P<:ParticleHilbertSpace} = speciescount(P)
# speciescount(::Type{ParticleHilbertSpace{PS, BR, QN}}) where {PS, BR, QN} = speciescount(PS)
# getspecies(::Type{ParticleHilbertSpace{PS, BR, QN}}, args...) where {PS, BR, QN} = getspecies(PS, args...)


for fname in [
    :exchangesign,
    :numspecies, :speciescount, :getspecies, :getspeciesname,
]
    @eval begin
        ($fname)(::Type{ParticleHilbertSpace{PS, BR, QN}}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
        ($fname)(::ParticleHilbertSpace{PS, BR, QN}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
    end
end




function bitoffset(phs::ParticleHilbertSpace{PS, BR, QN}, iptl::Integer, isite::Integer) where {PS, BR, QN}
    return phs.bitoffsets[isite] + bitoffset(PS, iptl)
end

function bitoffset(phs::ParticleHilbertSpace{PS, BR, QN}, isite::Integer) where {PS, BR, QN}
    return phs.bitoffsets[isite]
end


function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    iptl::Integer,
    isite::Integer,
)::BR where {PS, BR, QN}
    @boundscheck !(1<= iptl <= speciescount(PS)) && throw(BoundsError(PS.parameters, iptl))
    @boundscheck !(1<= isite <= length(phs.sites)) && throw(BoundsError(phs.sites, isite))
    n_sites = length(phs.sites)
    bm = get_bitmask(PS, iptl, BR)
    return bm << phs.bitoffsets[isite]
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    iptls::AbstractVector{<:Integer},
    isites::AbstractVector{<:Integer},
)::BR where {PS, BR, QN}
    @boundscheck for iptl in iptls
        !(1<= iptl <= speciescount(PS)) && throw(BoundsError(PS.parameters, iptl))
    end
    @boundscheck for isite in isites
        !(1<= isite <= length(phs.sites)) && throw(BoundsError(phs.sites, isite))
    end
    bm = mapreduce(iptl -> get_bitmask(PS, iptl, BR), |, iptls; init=zero(BR))
    return mapreduce(isite -> (bm << phs.bitoffsets[isite]), |, isites; init=zero(BR))
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    iptl::Integer,
    isites::AbstractVector{<:Integer},
)::BR where {PS, BR, QN}
    @boundscheck !(1<= iptl <= speciescount(PS)) && throw(BoundsError(PS.parameters, iptl))
    @boundscheck for isite in isites
        !(1<= isite <= length(phs.sites)) && throw(BoundsError(phs.sites, isite))
    end
    bm = get_bitmask(PS, iptl, BR)
    return mapreduce(isite -> (bm << phs.bitoffsets[isite]), |, isites; init=zero(BR))
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    iptls::AbstractVector{<:Integer},
    isite::Integer,
)::BR where {PS, BR, QN}
    @boundscheck for iptl in iptls
        !(1<= iptl <= speciescount(PS)) && throw(BoundsError(PS.parameters, iptl))
    end
    @boundscheck !(1<= isite <= length(phs.sites)) && throw(BoundsError(phs.sites, isite))
    bm = mapreduce(iptl -> get_bitmask(PS, iptl, BR), |, iptls; init=zero(BR))
    return bm << phs.bitoffsets[isite]
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    iptl::Integer,
    ::Colon,
)::BR where {PS, BR, QN}
    @boundscheck !(1<= iptl <= speciescount(PS)) && throw(BoundsError(PS.parameters, iptl))
    n_sites = length(phs.sites)
    bm = get_bitmask(PS, iptl, BR)
    return mapreduce((isite) -> bm << phs.bitoffsets[isite], |, 1:n_sites; init=zero(BR))
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    ::Colon,
    isite::Integer,
)::BR where {PS, BR, QN}
    @boundscheck !(1<= isite <= length(phs.sites)) && throw(BoundsError(phs.sites, isite))
    return make_bitmask(phs.bitoffsets[isite+1], phs.bitoffsets[isite], BR)
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    ::Colon,
    ::Colon,
)::BR where {PS, BR, QN}
    return make_bitmask(bitwidth(phs), BR)
end

# colon & vector.


function get_bitmask(phs::ParticleHilbertSpace{PS, BR, QN})::BR where {PS, BR, QN}
    return make_bitmask(bitwidth(phs), BR)
end

function get_occupancy(
    hs::ParticleHilbertSpace,
    iptl::Integer,
    isite::Integer,
    bvec::Unsigned,
)
    bm = get_bitmask(hs, iptl, isite)
    return Int( (bm & bvec) >> bitoffset(hs, iptl, isite) )
end


function set_occupancy(
    hs::ParticleHilbertSpace{PS, BR, QN},
    iptl::Integer,
    isite::Integer,
    bvec::BR2,
    count::Integer,
) where {PS, BR, QN, BR2}
    @boundscheck !(0 <= count <= maxoccupancy(getspecies(PS, iptl))) && throw(ArgumentError("count out of bounds"))
    bm = get_bitmask(hs, iptl, isite)
    return (bvec & ~bm) | (BR2(count) << bitoffset(hs, iptl, isite))
end


function quantum_number_sectors(phs::ParticleHilbertSpace{PS, BR, QN})::Vector{QN} where {PS, BR, QN}
    qns = Set{QN}([tuplezero(QN)])
    for site in phs.sites
        qns_next = Set{QN}()
        for state in site.states, q in qns
            push!(qns_next, q .+ state.quantum_number)
        end
        qns = qns_next
    end
    return sort(collect(qns))
end


function get_quantum_number(hs::ParticleHilbertSpace, binrep::Unsigned)
    return mapreduce(
        identity,
        tupleadd,
        let b = (get_bitmask(hs, :, isite) & binrep) >> bitoffset(hs, isite)
            #i = get_state_index(site, binrep)
            #site.states[i].quantum_number
            get_state(site, b).quantum_number
        end
            for (isite, site) in enumerate(hs.sites)
    )
end


function get_quantum_number(hs::ParticleHilbertSpace, statevec::AbstractVector{<:Integer})
    return mapreduce(identity, tupleadd,
        site.states[statevec[isite]].quantum_number
        for (isite, site) in enumerate(hs.sites)
    )
end


#
# function extract(hs::ParticleHilbertSpace{PS, BR, QN}, binrep::Unsigned)::CartesianIndex where {PS, BR, QN}
#   out = Int[]
#   for (isite, site) in enumerate(hs.sites)
#     @inbounds mask = make_bitmask(hs.bitwidths[isite], BR)
#     index = Int(binrep & mask) + 1
#     @boundscheck if !(1 <= index <= length(site.states))
#       throw(BoundsError(1:length(site.states), index))
#     end
#     push!(out, index)
#     binrep = binrep >> hs.bitwidths[isite]
#   end
#   return CartesianIndex(out...)
# end







#=
export canonize
function canonize(pt::Type{<:Fermion}, locations::AbstractVector{<:Integer})
  sgn = 1
  out = Int.(locations)
  n = length(locations)
  for i in n:-1:1
    for j in Base.OneTo(i-1)
      lhs, rhs = out[j], out[j+1]
      if lhs > rhs
        sgn = -sgn
        out[j], out[j+1] = rhs, lhs
      end
    end
  end
  return (out, sgn)
end

function canonize(pt::Type{<:Boson}, locations::AbstractVector{<:Integer})
    return (sort(locations), 1)
end

function canonize(pt::Type{<:HardcoreBoson}, locations::AbstractVector{<:Integer})
    return (sort(locations), 1)
end



=#
