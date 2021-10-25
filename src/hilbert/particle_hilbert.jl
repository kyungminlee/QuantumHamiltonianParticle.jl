using QuantumHamiltonian
import QuantumHamiltonian.get_bitmask

export ParticleHilbertSpace

export get_parity_bitmask
export get_fermion_parity
export get_occupancy, set_occupancy

import QuantumHamiltonian.AbstractHilbertSpace

import QuantumHamiltonian.scalartype
import QuantumHamiltonian.qntype
import QuantumHamiltonian.basespace
import QuantumHamiltonian.numsites
import QuantumHamiltonian.get_site
import QuantumHamiltonian.bitwidth
import QuantumHamiltonian.bitoffset
import QuantumHamiltonian.get_bitmask
import QuantumHamiltonian.quantum_number_sectors
import QuantumHamiltonian.get_quantum_number
import QuantumHamiltonian.compress
import QuantumHamiltonian.extract
import QuantumHamiltonian.uncompress


# Add a decoration to the existing Hilbert space
"""
    ParticleHilbertSpace{PS, BR, QN} <: AbstractHilbertSpace{QN}

Particle Hilbert space.

# Example
```
ParticleHilbertSpace([site1, site2, site3, site4])
```
"""
struct ParticleHilbertSpace{
        PS<:ParticleSector,
        BR<:Unsigned,
        QN<:Tuple{Vararg{<:AbstractQuantumNumber}}
}<:AbstractHilbertSpace{QN}
    sites::Vector{ParticleSite{PS, BR, QN}}

    bitwidths::Matrix{Int}
    bitoffsets::Matrix{Int}
    bitmasks::Matrix{BR}
    parity_bitmasks::Matrix{BR}

    function ParticleHilbertSpace(sites::AbstractVector{ParticleSite{PS, BR, QN}}) where {PS, BR, QN}
        nsites = length(sites)
        nptls = numspecies(PS)

        bitwidths = Matrix{Int}(undef, (nptls+1, nsites))
        bitoffsets = Matrix{Int}(undef, (nptls+1, nsites+1))
        bitmasks = Matrix{BR}(undef, (nptls+1, nsites+1))

        offset = 0
        for isite in 1:nsites
            for iptl in 1:nptls
                bw = bitwidth(PS, iptl)
                bitwidths[iptl, isite] = bw
                bitoffsets[iptl, isite] = offset

                bitmasks[iptl, isite] = make_bitmask(offset + bw, offset, BR)
                offset += bw
            end
            sbw = sum(view(bitwidths, 1:nptls, isite))
            bitwidths[end, isite] = sbw
            bitoffsets[end, isite] = bitoffsets[1, isite]
            bitmasks[end, isite] = make_bitmask(bitwidths[end, isite] + bitoffsets[1, isite], bitoffsets[1, isite], BR)
        end
        bitoffsets[:, end] .= offset
        for iptl in 1:nptls
            bitmasks[iptl, end] = mapreduce(identity, |, view(bitmasks, iptl, 1:nsites); init=zero(BR))
        end
        bitmasks[end, end] = make_bitmask(bitoffsets[end, end], BR)
        if sizeof(BR) * 8 < bitoffsets[end, end]
            throw(ArgumentError("type $BR too small to represent the hilbert space (need $(bitoffsets[end]) bits)"))
        end
        parity_bitmasks = zeros(BR, (nptls, nsites))
        for iptl in 1:nptls
            species = getspecies(PS, iptl)
            if isboson(species) || isspin(species)
                # do nothing
            elseif isfermion(species)
                for isite in 1:nsites
                    bm_species = bitmasks[iptl, end]
                    bm_site = bitmasks[iptl, isite]
                    bm_mask = (bm_site - 1) & bm_species  # σᶻ in jordan wigner string
                    parity_bitmasks[iptl, isite] = bm_mask
                end
            else
                throw(ArgumentError("Particle species $(species) currently unsupported"))
            end
        end
        return new{PS, BR, QN}(sites, bitwidths, bitoffsets, bitmasks, parity_bitmasks)
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

numsites(hs::ParticleHilbertSpace) = length(hs.sites)

get_site(hs::ParticleHilbertSpace, i::Integer) = hs.sites[i]


for fname in [
    :exchangesign,
    :numspecies, :speciescount, :getspecies, :getspeciesname,
]
    @eval begin
        ($fname)(::Type{ParticleHilbertSpace{PS, BR, QN}}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
        ($fname)(::ParticleHilbertSpace{PS, BR, QN}, args...) where {PS, BR, QN} = ($fname)(PS, args...)
    end
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


"""
    get_quantum_number(phs, bitrep)

Get quantum number of the basis state `bitrep`.
"""
function get_quantum_number(hs::ParticleHilbertSpace{PS, BR, QN}, binrep::Unsigned) where {PS, QN, BR}
    bw = bitwidth(PS)
    return mapreduce(
        identity,
        tupleadd,
        let b = (get_bitmask(hs, :, isite) & binrep) >> (bw*(isite-1))
            get_state(site, b).quantum_number
        end
            for (isite, site) in enumerate(hs.sites)
    )
end


"""
    get_quantum_number(phs, statevec)

Get quantum number of the basis state `statevec`.
"""
function get_quantum_number(hs::ParticleHilbertSpace, statevec::AbstractVector{<:Integer})
    return mapreduce(identity, tupleadd,
        site.states[statevec[isite]].quantum_number
        for (isite, site) in enumerate(hs.sites)
    )
end


"""
    bitwidth(phs)

Return number of bits needed to represent basis states of `phs`.
"""
bitwidth(hs::ParticleHilbertSpace) = hs.bitoffsets[end, end]

function bitwidth(phs::ParticleHilbertSpace, isite::Integer)
    nsites = numsites(phs)
    @boundscheck !(1 <= isite <= nsites) && throw(BoundsError(phs.sites, isite))
    @inbounds return phs.bitwidths[end, isite]
end

function bitwidth(phs::ParticleHilbertSpace, iptl::Integer, isite::Integer)
    nptls = numspecies(phs)
    nsites = numsites(phs)
    @boundscheck !(1 <= iptl <= nptls) && throw(BoundsError(PS.parameters, iptl))
    @boundscheck !(1 <= isite <= nsites) && throw(BoundsError(phs.sites, isite))
    @inbounds return phs.bitwidths[iptl, isite]
end

bitoffset(hs::ParticleHilbertSpace, isite::Integer) = hs.bitoffsets[1, isite]

"""
    bitoffset(phs, iptl, isite)

Get the bit offset of the particle `iptl` at site `isite`.
"""
function bitoffset(phs::ParticleHilbertSpace{PS, BR, QN}, iptl::Integer, isite::Integer) where {PS, BR, QN}
    nptls = numspecies(phs)
    nsites = numsites(phs)
    @boundscheck !(1 <= iptl <= nptls) && throw(BoundsError(PS.parameters, iptl))
    @boundscheck !(1 <= isite <= nsites) && throw(BoundsError(phs.sites, isite))
    @inbounds return phs.bitoffsets[iptl, isite]
end





"""
    get_bitmask(phs, iptl, isite)

Get the bit mask for the particles `iptl` at sites `isite`.
`iptl` or `isite` can either be integer, a vector of integers, or colon (:).
Bitwise or is taken over list of iptl.
"""
function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    iptl::Integer,
    isite::Integer,
)::BR where {PS, BR, QN}
    nptls = numspecies(phs)
    nsites = numsites(phs)
    @boundscheck !(1 <= iptl <= nptls) && throw(BoundsError(PS.parameters, iptl))
    @boundscheck !(1 <= isite <= nsites) && throw(BoundsError(phs.sites, isite))
    @inbounds return phs.bitmasks[iptl, isite]
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    iptls::AbstractVector{<:Integer},
    isites::AbstractVector{<:Integer},
)::BR where {PS, BR, QN}
    nptls = numspecies(phs)
    nsites = numsites(phs)
    @boundscheck for iptl in iptls
        !(1 <= iptl <= nptls) && throw(BoundsError(PS.parameters, iptl))
    end
    @boundscheck for isite in isites
        !(1 <= isite <= nsites) && throw(BoundsError(phs.sites, isite))
    end
    bm = zero(BR)
    for iptl in iptls, isite in isites
        @inbounds bm |= phs.bitmasks[iptl, isite]
    end
    return bm
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    iptl::Integer,
    isites::AbstractVector{<:Integer},
)::BR where {PS, BR, QN}
    nptls = numspecies(phs)
    nsites = numsites(phs)
    @boundscheck !(1 <= iptl <= nptls) && throw(BoundsError(PS.parameters, iptl))
    @boundscheck for isite in isites
        !(1 <= isite <= nsites) && throw(BoundsError(phs.sites, isite))
    end
    # bm = get_bitmask(PS, iptl, BR)
    # return mapreduce(isite -> (bm << phs.bitoffsets[isite]), |, isites; init=zero(BR))
    bm = zero(BR)
    for isite in isites
        @inbounds bm |= phs.bitmasks[iptl, isite]
    end
    return bm
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    iptls::AbstractVector{<:Integer},
    isite::Integer,
)::BR where {PS, BR, QN}
    nptls = numspecies(phs)
    nsites = numsites(phs)
    @boundscheck for iptl in iptls
        !(1 <= iptl <= nptls) && throw(BoundsError(PS.parameters, iptl))
    end
    @boundscheck !(1 <= isite <= nsites) && throw(BoundsError(phs.sites, isite))
    bm = zero(BR)
    for iptl in iptls
        @inbounds bm |= phs.bitmasks[iptl, isite]
    end
    return bm    
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    iptl::Integer,
    ::Colon,
)::BR where {PS, BR, QN}
    nptls = numspecies(phs)
    @boundscheck !(1 <= iptl <= nptls) && throw(BoundsError(PS.parameters, iptl))
    @inbounds return phs.bitmasks[iptl, end]
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    ::Colon,
    isite::Integer,
)::BR where {PS<:ParticleSector, BR, QN}
    nsites = numsites(phs)
    @boundscheck !(1 <= isite <= nsites) && throw(BoundsError(phs.sites, isite))
    @inbounds return phs.bitmasks[end, isite]
end

function get_bitmask(
    phs::ParticleHilbertSpace{PS, BR, QN},
    ::Colon,
    ::Colon,
)::BR where {PS, BR, QN}
    return phs.bitmasks[end, end]
end

function get_bitmask(phs::ParticleHilbertSpace{PS, BR, QN})::BR where {PS, BR, QN}
    return phs.bitmasks[end, end]
end


"""
    get_parity_bitmask(phs, iptl, isite)

Get parity bitmask (i.e. position of the Wigner-Jordan string). Nonzero only for fermions.
"""
function get_parity_bitmask(phs::ParticleHilbertSpace{PS, BR, QN}, iptl::Integer, isite::Integer) where {PS, BR, QN}
    nptls = numspecies(phs)
    nsites = numsites(phs)
    @boundscheck !(1 <= iptl <= nptls) && throw(BoundsError(PS.parameters, iptl))
    @boundscheck !(1 <= isite <= nsites) && throw(BoundsError(phs.sites, isite))
    @inbounds return phs.parity_bitmasks[iptl, isite]
end


"""
    get_occupancy(phs, iptl, isite, bvec::Unsigned)

Get occupancy of particle `iptl` at site `isite` for the given basis state `bvec`.
"""
function get_occupancy(
    hs::ParticleHilbertSpace,
    iptl::Integer,
    isite::Integer,
    bvec::Unsigned,
)
    bm = get_bitmask(hs, iptl, isite)
    return Int( (bm & bvec) >> bitoffset(hs, iptl, isite) )
end


"""
    set_occupancy(phs, iptl, isite, bvec::Unsigned, count)

Set occupancy of particle `iptl` at site `isite` for the given basis state `bvec` to `count`.
"""
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


"""
    compress(hs, indexarray, [type])

Return the binary representation of the basis state represented by `indexarray`, optionally in type `type`.

# Arguments
- `hs::ParticleHilbertSpace{PS, BR, QN}`
- `indexarray::CartesianIndex`
- `type::Type{BR2}=BR`
"""
function compress(
    hs::ParticleHilbertSpace{PS, BR, QN},
    indexarray::CartesianIndex,
    ::Type{BR2}=BR,
) where {PS, BR, QN, BR2<:Unsigned}
    return BR2(statevec2occbin(hs, collect(indexarray.I)))
end


"""
    extract(hs, occbin)

Return the CartesianIndex representation of the basis state represented by `occbin`.

# Arguments
- `hs::ParticleHilbertSpace`
- `occbin::Unsigned`
"""
function extract(hs::ParticleHilbertSpace, occbin::Unsigned)
    return CartesianIndex(occbin2statevec(hs, occbin)...)
end


"""
    uncompress(hs, occbin)

Return the CartesianIndex representation of the basis state represented by `occbin`.

# Arguments
- `hs::ParticleHilbertSpace`
- `occbin::Unsigned`
"""
function uncompress(hs::ParticleHilbertSpace, occbin::Unsigned)
    return CartesianIndex(occbin2statevec(hs, occbin)...)
end


function Base.keys(hs::ParticleHilbertSpace)
    return CartesianIndices(((1:length(site.states) for site in hs.sites)...,))
end


function Base.getindex(
    hs::ParticleHilbertSpace{PS, BR, QN},
    idx::CartesianIndex,
) where {PS, BR, QN}
    return compress(hs, idx, BR)
end


function Base.getindex(
    hs::ParticleHilbertSpace{PS, BR, QN},
    idx::Integer...
) where {PS, BR, QN}
    return compress(hs, CartesianIndex(idx...), BR)
end
