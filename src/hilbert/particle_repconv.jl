# conversion between representations

# Single Sites (ParticleSite)
# ---------------------------
# OccBin <-> OccVec <-> State
#
# Multiple Sites (ParticleHilbertSpace)
# -------------------------------------
# OccBin <-> OccMat <-> StateVec <-> LocVec
#                                <-(has sign)
#                                -> (no sign)

# ------------------
#   ParticleSite
# ------------------


export occbin2state
export occbin2occvec

export occvec2occbin
export occvec2state

export state2occbin
export state2occvec

function occbin2occvec(ps::ParticleSite{PS, BR, QN}, occbin::Unsigned) where {PS, BR, QN}
    return extract(PS, occbin)
end

function occvec2occbin(ps::ParticleSite{PS, BR, QN}, occarr::AbstractVector{<:Integer}) where {PS, BR, QN}
    return compress(PS, occarr, BR)
end

function state2occbin(ps::ParticleSite{PS, BR, QN}, index::Integer)::BR where {PS, BR, QN}
    return ps.states[index].occupancy_binary
end

function state2occvec(ps::ParticleSite{PS, BR, QN}, index::Integer)::Vector{Int} where {PS, BR, QN}
    return extract(PS, ps.states[index].occupancy_binary)
end

function occbin2state(ps::ParticleSite{PS, BR, QN}, occbin::Unsigned)::Int where {PS, BR, QN}
    return ps.state_lookup[occbin]
end

function occvec2state(ps::ParticleSite{PS, BR, QN}, occarr::AbstractVector{<:Integer})::Int where {PS, BR, QN}
    return ps.state_lookup[compress(PS, occarr)]
end


# --------------------
# ParticleHilbertSpace
# --------------------
export statevec2occbin
export statevec2occmat
export statevec2locvec

export occmat2occbin
export occmat2statevec
export occmat2locvec

export occbin2occmat
export occbin2statevec
export occbin2locvec

export locvec2occmat
export locvec2occbin
export locvec2statevec


function statevec2occbin(phs::ParticleHilbertSpace, statevec::CartesianIndex)
    return statevec2occbin(phs, collect(statevec.I))
end


function statevec2occmat(phs::ParticleHilbertSpace, statevec::CartesianIndex)
    return statevec2occmat(phs, collect(statevec.I))
end


function statevec2locvec(phs::ParticleHilbertSpace, statevec::CartesianIndex)
    return statevec2locvec(phs, collect(statevec.I))
end


function statevec2occbin(
    phs::ParticleHilbertSpace{PS, BR, QN},
    statevec::AbstractVector{<:Integer},
)::BR where {PS, BR, QN}
    if length(statevec) != length(phs.sites)
        throw(ArgumentError("Length of statevec does not match the number of sites"))
    end
    occbin = zero(BR)
    for (isite, (site, istate)) in enumerate(zip(phs.sites, statevec))
        siteoccbin = state2occbin(site, istate)
        occbin |= siteoccbin << bitoffset(phs, isite)
    end
    return occbin
end


function statevec2occmat(
    ps::ParticleHilbertSpace{PS, BR, QN},
    statevec::AbstractVector{<:Integer},
)::Matrix{Int} where {PS, BR, QN}
    nptls = speciescount(PS)
    nsites = length(ps.sites)
    if length(statevec) != nsites
        throw(ArgumentError("Length of statevec does not match the number of sites"))
    end
    occarr = zeros(Int, (nptls, nsites))
    for (isite, (site, istate)) in enumerate(zip(ps.sites, statevec))
        occarr[:, isite] = extract(PS, state2occbin(site, istate))
    end
    return occarr
end


function statevec2locvec(
    phs::ParticleHilbertSpace{PS, BR, QN},
    statevec::AbstractVector{<:Integer},
) where {PS, BR, QN}
    nsites = length(phs.sites)
    nptls = speciescount(PS)
    occmat = statevec2occmat(phs, statevec)
    out = [Int[] for _ in 1:nptls]
    for iptl in 1:nptls, isite in 1:nsites
        append!(out[iptl], isite for c in 1:occmat[iptl, isite])
    end
    return out
end


function occmat2occbin(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occmat::AbstractMatrix{<:Integer},
) where {PS, BR, QN}
    nptls = speciescount(PS)
    nsites = length(phs.sites)
    size(occmat) != (nptls, nsites) && throw(ArgumentError("Wrong occmat size"))
    occbin = zero(BR)
    for isite in 1:nsites
        occbin |= compress(PS, view(occmat, :, isite), BR) << bitoffset(phs, isite)
    end
    return occbin
end


function occmat2statevec(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occmat::AbstractMatrix{<:Integer},
) where {PS, BR, QN}
    nptls = speciescount(PS)
    nsites = length(phs.sites)
    size(occmat) != (nptls, nsites) && throw(ArgumentError("Wrong occmat size"))
    out = Vector{Int}(undef, length(phs.sites))
    for (isite, site) in enumerate(phs.sites)
        siteoccbin = compress(PS, occmat[:,isite], BR)
        idx = occbin2state(site, siteoccbin)
        out[isite] = idx
    end
    return out
end


function occmat2locvec(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occmat::AbstractMatrix{<:Integer},
) where {PS, BR, QN}
    nsites = length(phs.sites)
    nptls = speciescount(PS)
    size(occmat) != (nptls, nsites) && throw(ArgumentError("Wrong occmat size"))
    out = [Int[] for i in 1:nptls]
    for iptl in 1:nptls, isite in 1:nsites
        append!(out[iptl], isite for c in 1:occmat[iptl, isite])
    end
    return out
end


function occbin2statevec(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occbin::Unsigned,
)::Vector{Int} where {PS, BR, QN}
    out = Vector{Int}(undef, length(phs.sites))
    for (isite, site) in enumerate(phs.sites)
        siteoccbin = (occbin & get_bitmask(phs, :, isite)) >> bitoffset(phs, isite)
        idx = occbin2state(site, siteoccbin)
        out[isite] = idx
    end
    return out
end


function occbin2occmat(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occbin::Unsigned,
) where {PS, BR, QN}
    nptls = speciescount(PS)
    nsites = length(phs.sites)
    occmat = Matrix{Int}(undef, (nptls, nsites))
    for isite in 1:nsites
        siteoccbin = (occbin & get_bitmask(phs, :, isite)) >> bitoffset(phs, isite)
        extract(PS, siteoccbin, view(occmat, :, isite))
    end
    return occmat
end


function occbin2occmat(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occbin::Unsigned,
    occmat::AbstractMatrix{<:Integer},
) where {PS, BR, QN}
    nsites = length(phs.sites)
    for isite in 1:nsites
        siteoccbin = (occbin & get_bitmask(phs, :, isite)) >> bitoffset(phs, isite)
        extract(PS, siteoccbin, view(occmat, :, isite))
    end
    return occmat
end


function occbin2locvec_naive(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occbin::Unsigned,
) where {PS, BR, QN}
    occmat = occbin2occmat(phs, occbin)
    return occmat2locvec(phs, occmat)
end


function occbin2locvec(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occbin::Unsigned,
) where {PS, BR, QN}
    nsites = length(phs.sites)
    nptls = speciescount(PS)
    out = [Int[] for _ in 1:nptls]
    for iptl in 1:nptls
        sizehint!(out[iptl], nsites)
    end
    sbw::Int = bitwidth(PS)
    for iptl in 1:nptls
        pbw::Int = bitwidth(PS, iptl)
        pbm::BR = make_bitmask(pbw, BR)
        for isite in 1:nsites
            occ = Int( (occbin >> (sbw * (isite-1))) & pbm )
            append!(out[iptl], isite for _ in 1:occ)
        end
        occbin >>= pbw
    end
    return out
end


function occbin2locvec!(
    out::AbstractVector{<:AbstractVector{<:Integer}},
    phs::ParticleHilbertSpace{PS, BR, QN},
    occbin::Unsigned,
) where {PS, BR, QN}
    nsites = length(phs.sites)
    nptls = speciescount(PS)
    sbw::Int = bitwidth(PS)
    for iptl in 1:nptls
        pbw::Int = bitwidth(PS, iptl)
        pbm::BR = make_bitmask(pbw, BR)
        for isite in 1:nsites
            occ = Int( (occbin >> (sbw * (isite-1))) & pbm )
            append!(out[iptl], isite for _ in 1:occ)
        end
        occbin >>= pbw
    end
    return out
end

function locvec2occmat(
    phs::ParticleHilbertSpace{PS, BR, QN},
    particle_locations::AbstractVector{<:AbstractVector{<:Integer}},
) where {PS, BR, QN}
    nsites = length(phs.sites)
    nptls = speciescount(PS)
    if length(particle_locations) != nptls
        throw(ArgumentError("length of particle locations does not match the number of particles"))
    end
    sgn = 1
    occmat = zeros(Int, (nptls, nsites)) # occupation matrix
    for (iptl, p) in enumerate(particle_locations)
        if exchangesign(PS, iptl) == 1
            # do nothing
        elseif exchangesign(PS, iptl) == -1
            sgn = isparityodd(p) ? -sgn : sgn
        else
            throw(ArgumentError("Unsupported exchange sign $(exchangesign(getspecies(PS, iptl))) for $(getspecies(PS, iptl))")) # COV_EXCL_LINE
        end
        for isite in p
            occmat[iptl, isite] += 1
        end
        # should not do `occmat[iptl, p] .+= 1`  (may contain duplicates)
    end
    return (occmat, sgn)
end



function locvec2occbin_naive(
    phs::ParticleHilbertSpace{PS, BR, QN},
    particle_locations::AbstractVector{<:AbstractVector{<:Integer}},
) where {PS, BR, QN}
    occmat, sgn = locvec2occmat(phs, particle_locations)
    return (occmat2occbin(phs, occmat), sgn)
end


function locvec2occbin(
    phs::ParticleHilbertSpace{PS, BR, QN},
    particle_locations::AbstractVector{<:AbstractVector{<:Integer}},
) where {PS, BR, QN}
    nptls = speciescount(PS)
    if length(particle_locations) != nptls
        throw(ArgumentError("length of particle locations does not match the number of particles"))
    end
    sgn = 1
    occbin = zero(BR)
    for (iptl, p) in enumerate(particle_locations)
        if exchangesign(PS, iptl) == 1
            # do nothing
        elseif exchangesign(PS, iptl) == -1
            sgn = isparityodd(p) ? -sgn : sgn
        else
            throw(ArgumentError("Unsupported exchange sign $(exchangesign(getspecies(PS, iptl))) for $(getspecies(PS, iptl))")) # COV_EXCL_LINE
        end
        for isite in p
            count = get_occupancy(PS, iptl, isite, occbin)
            occbin = set_occupancy(PS, iptl, isite, occbin, count+1)
        end
    end
    return (occbin, sgn)
end


"""
    locvec2statevec


particles : Vector of (Vector of particle location)

# Example

```
  P S 1 2 3 4 5 = |0, e↓, e↑, m, m⟩
  e↑  0 0 1 0 0
  e↓  0 1 0 0 0 = c†(m,4) c†(m,5) c†(e↑,3) c†(e↓,2) |Ω⟩
  m   0 0 0 1 1
```
"""
function locvec2statevec(
    phs::ParticleHilbertSpace{PS, BR, QN},
    particle_locations::AbstractVector{<:AbstractVector{<:Integer}},
) where {PS, BR, QN}
    occmat, sgn = locvec2occmat(phs, particle_locations)
    return (occmat2statevec(phs, occmat), sgn)
end
