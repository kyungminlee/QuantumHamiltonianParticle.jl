export ParticleHilbertSpace
export getparticlebitstring
export getparticlelist
export getparticlelocationlist

import ExactDiagonalization.extract


struct ParticleHilbertSpace{
    P<:ParticleSector,
    QN<:Tuple{Vararg{<:AbstractQuantumNumber}},
}<:AbstractHilbertSpace

    sites::Vector{ParticleSite{P, QN}}
    bitwidths::Vector{Int}
    bitoffsets::Vector{Int}

    function ParticleHilbertSpace(sites::AbstractVector{ParticleSite{P, QN}}) where {P, QN}
        bitwidths = map(bitwidth, sites)
        bitoffsets = Int[0, cumsum(bitwidths)...]
        new{P, QN}(sites, bitwidths, bitoffsets)
    end
end


function extract(
    hs::ParticleHilbertSpace{P, QN},
    binrep::BR
)::CartesianIndex where {P, QN, BR<:Unsigned}
    out = Int[]
    for (isite, site) in enumerate(hs.sites)
        @inbounds mask = make_bitmask(hs.bitwidths[isite], BR)
        index = Int(binrep & mask) + 1
        @boundscheck if !(1 <= index <= length(site.states))
            throw(BoundsError(1:length(site.states), index))
        end
        push!(out, index)
        binrep = binrep >> hs.bitwidths[isite]
    end
    return CartesianIndex(out...)
end


function getparticlelist(hs::ParticleHilbertSpace, idx::CartesianIndex)
    return [site.states[i].particles for (site, i) in zip(hs.sites, idx.I)]
end


function getparticlelocationlist(
    hs::ParticleHilbertSpace{P, QN},
    idx::CartesianIndex,
) where {P, QN}
    out = [Int[] for i in 1:numparticletype(P)]  # TODO: cleanup with a method
    for (isite, (site, istate)) in enumerate(zip(hs.sites, idx.I))
        state = site.states[istate]
        for (iparticle, particle) in enumerate(state.particlesector.particlecounts)
            append!(out[iparticle], isite for _ in 1:particle.value)
        end
    end
    return out
end


export getparticleoccupation

function getparticleoccupation(
    phs::ParticleHilbertSpace{P, QN},
    idx::CartesianIndex,
) where {P, QN}
    out = Matrix{Int}(undef, (particletypecount(P), length(phs.sites)))
    for (isite, (site, istate)) in enumerate(zip(phs.sites, idx.I))
        state = site.states[istate]
        for (iparticle, particle) in enumerate(state.particlesector.particlecounts)
            out[iparticle, isite] = particle.value
        end
    end
    return out
end


function getparticlebitstring(
    hs::ParticleHilbertSpace{ParticleSector{P}, QN},
    idx::CartesianIndex,
    ::Type{BR},
) where {P, QN, BR<:Unsigned}
    bitwidths = particlebitwidth.(P.parameters)
    bitshifts = [0, cumsum(bitwidths)...]
    bitwidth = sum(bitwidths)
    out = zero(BR)
    for (isite, (site, istate)) in enumerate(zip(hs.sites, idx.I))
        state = site.states[istate]
        for (particlecount, shift) in zip(state.particlesector.particlecounts, bitshifts)
            out |= BR(particlecount.value) << ((isite-1)*bitwidth + shift)
        end
    end
    return out
end


function getparticlebitstring(
    hs::ParticleHilbertSpace{ParticleSector{P}, QN},
    bvec::BR1,
    ::Type{BR},
) where {P, QN, BR1<:Unsigned, BR<:Unsigned}
    return getparticlebitstring(hs, extract(hs, bvec), BR)
end
