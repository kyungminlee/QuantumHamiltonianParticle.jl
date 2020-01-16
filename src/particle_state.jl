export State, Site, HilbertSpace
export particle_decompose
export particle_compose

export StateParticleMap

using ExactDiagonalization


struct StateParticleMap{NP} # NP: number of particle types
    state_to_particle::Matrix{Int}
    particle_to_state::Dict{NTuple{NP, Int}, Int}
    function StateParticleMap(stateocc::AbstractVector{T}) where {T<:NTuple{NP, <:Integer}} where {NP}
        n_states = length(stateocc)
        s2p = Matrix{Int}(undef, (NP, n_states))
        p2s = Dict{NTuple{NP, Int}, Int}()
        for (i, occ) in enumerate(stateocc)
            @assert !haskey(p2s, occ)
            # check all positive
            if any(x < 0 for x in occ)
                throw(ArgumentError("Particle occupation need to be positive"))
            end
            s2p[:, i] = collect(occ)
            p2s[tuple(occ...)] = i
        end
        return new{NP}(s2p, p2s)
    end
end



export getspeciessymb
export maxoccupancy
export exchangesign

getspeciessymb(p::Type{Boson{S, M}}) where {S, M} = S
getspeciessymb(p::Type{HardcoreBoson{S}}) where S = S
getspeciessymb(p::Type{Fermion{S}}) where S = S

# Default is 1 for different types
exchangesign(::Type{<:AbstractParticle}) = 1
exchangesign(::Type{<:Fermion}) = -1
exchangesign(::T) where {T<:AbstractParticle} = exchangesign(T)

maxoccupancy(::Type{<:Fermion}) = 1
maxoccupancy(::Type{<:HardcoreBoson}) = 1
maxoccupancy(::Type{Boson{S, M}}) where {S, M} = M::Int

exchangesign(pt::Type{<:Fermion}, locations::AbstractVector{<:Integer}) = parity(locations) ? -1 : 1
exchangesign(pt::Type{<:Boson}, locations::AbstractVector{<:Integer}) = 1
exchangesign(pt::Type{<:HardcoreBoson}, locations::AbstractVector{<:Integer}) = 1

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


# in utility
function parity(arr::AbstractVector{T}) where T
    n = length(arr)
    parity = true
    for i in n-1:-1:1
        x = arr[i]
        for j in i+1:n
            if arr[j] < x
                parity = !parity
            end
        end
    end
    return parity
end



export ParticleHilbertSpace

# Add a decoration to the existing Hilbert space
struct ParticleHilbertSpace{HS, NP} #<:AbstractHilbertSpace
    parent::HS
    particle_types::NTuple{NP, DataType}
    state_particle_maps::Vector{StateParticleMap{NP}}

    function ParticleHilbertSpace(hs::HS,
        particle_types,
        spms::AbstractVector{StateParticleMap{NP}}
        ) where {HS, NP}
        
        #TODO check a couple of things
        particle_types_tuple = tuple(particle_types...)
        return new{HS, NP}(hs, particle_types_tuple, spms)
    end
end



export particle_decompose
"""
    particle_decompose


"""
function particle_decompose(
        ph::ParticleHilbertSpace{HS, NP},
        indexarray::AbstractVector{<:Integer}) where {HS, NP}

    maps = ph.state_particle_maps

    if length(maps) != length(indexarray)
        throw(ArgumentError("map and indexarray should have the same length"))
    end

    out = [Int[] for i in Base.OneTo(NP)]
    for (isite, (map, istate)) in enumerate(zip(maps, indexarray))
        for iptl in 1:NP
            append!(out[iptl], isite for c in 1:map.state_to_particle[iptl, istate])
        end
    end
    return out
end


export particle_compose

"""
    particle_compose


particles : Vector of (Vector of particle location)

# Example

```
  P S 1 2 3 4 5 = |0, e↓, e↑, m, m⟩
  e↑  0 0 1 0 0
  e↓  0 1 0 0 0 = c†(m,4) c†(m,5) c†(e↑,3) c†(e↓,2) |Ω⟩
  m   0 0 0 1 1
```
"""
function particle_compose(
        ph::ParticleHilbertSpace{HS, NP},
        particles::AbstractVector{<:AbstractVector{<:Integer}}) where {HS, NP}

    maps = ph.state_particle_maps
    nsites = length(maps)
    occmat = zeros(Int, (NP, nsites)) # occupation matrix
    sgn = 1
    for (iptl, p) in enumerate(particles)
        sgn *= exchangesign(ph.particle_types[iptl], p)
        for isite in p
            occmat[iptl, isite] += 1
        end
    end

    return (Int[
        let occtuple = tuple(occmat[:, isite]...) # occupation number tuple
            maps[isite].particle_to_state[ occtuple ]
        end
        for isite in 1:nsites], sgn)
end


# --- now the operators
















