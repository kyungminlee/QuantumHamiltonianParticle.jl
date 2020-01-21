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


export occupancy_binary_to_state
export occupancy_binary_to_occupancy_array

export occupancy_array_to_occupancy_binary
export occupancy_array_to_state

export state_to_occupancy_binary
export state_to_occupancy_array

occupancy_binary_to_occupancy_array = extract
occupancy_array_to_occupancy_binary = compress

function state_to_occupancy_binary(ps::ParticleSite{PS, BR, QN}, index::Integer)::BR where {PS, BR, QN}
  return ps.states[index].occupancy_binary
end

function state_to_occupancy_array(ps::ParticleSite{PS, BR, QN}, index::Integer)::Vector{Int} where {PS, BR, QN}
  return extract(PS, ps.states[index].occupancy_binary)
end

function occupancy_binary_to_state(ps::ParticleSite{PS, BR, QN}, occbin::Unsigned)::Int where {PS, BR, QN}
  return ps.state_lookup[occbin]
end

function occupancy_array_to_state(ps::ParticleSite{PS, BR, QN},
                                  occarr::AbstractVector{<:Integer})::Int where {PS, BR, QN}
  return ps.state_lookup[compress(PS, occarr)]
end


# --------------------
# ParticleHilbertSpace
# --------------------
export state_array_to_occupancy_binary
export state_array_to_occupancy_array
export state_array_to_location_array

export occupancy_array_to_occupancy_binary
export occupancy_array_to_state_array

export occupancy_binary_to_occupancy_array
export occupancy_binary_to_state_array

export location_array_to_state_array


## new
export occupancy_binary_to_location_array
export occupancy_array_to_location_array


export location_array_to_occupancy_array
export location_array_to_occupancy_binary



function state_array_to_occupancy_binary(
    phs::ParticleHilbertSpace{PS, BR, QN},
    indexarray::AbstractVector{<:Integer})::BR where {PS, BR, QN}
  if length(indexarray) != length(phs.sites)
    throw(ArgumentError("Length of indexarray does not match the number of sites"))
  end
  occbin = zero(BR)
  for (isite, (site, istate)) in enumerate(zip(phs.sites, indexarray))
    siteoccbin = state_to_occupancy_binary(site, istate)
    occbin |= siteoccbin << phs.bitoffsets[isite]
  end
  return occbin
end


function state_array_to_occupancy_array(
    ps::ParticleHilbertSpace{PS, BR, QN},
    indexarray::AbstractVector{<:Integer})::Matrix{Int} where {PS, BR, QN}

  n_particles = num_particle_species(PS)
  n_sites = length(ps.sites)

  if length(indexarray) != length(ps.sites)
    throw(ArgumentError("Length of indexarray does not match the number of sites"))
  end

  occarr = zeros(Int, (n_particles, n_sites))
  for (isite, (site, istate)) in enumerate(zip(ps.sites, indexarray))
    occarr[:, isite] = extract(PS, state_to_occupancy_binary(site, istate))
  end
  return occarr
end



function state_array_to_location_array(
    phs::ParticleHilbertSpace{PS, BR, QN},
    indexarray::AbstractVector{<:Integer}) where {PS, BR, QN}
  n_sites = length(phs.sites)
  n_particles = num_particle_species(PS)
  occmat = state_array_to_occupancy_array(phs, indexarray)
  out = [Int[] for i in 1:n_particles]
  for iptl in 1:n_particles
    for isite in 1:n_sites
      append!(out[iptl], isite for c in 1:occmat[iptl, isite])
    end
  end
  return out
end




function occupancy_array_to_occupancy_binary(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occmat::Matrix{<:Integer}) where {PS, BR, QN}
  n_particles = num_particle_species(PS)
  n_sites = length(phs.sites)
  if size(occmat) != (n_particles, n_sites)
    throw(ArgumentError("Wrong occmat size"))
  end
  occbin = zero(BR)
  for isite in 1:n_sites
    occbin |= compress(PS, occmat[:, isite]) << phs.bitoffsets[isite]
  end
  return occbin
end



function occupancy_array_to_state_array(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occmat::Matrix{<:Integer}) where {PS, BR, QN}
  n_particles = num_particle_species(PS)
  n_sites = length(phs.sites)

  if size(occmat) != (n_particles, n_sites)
    throw(ArgumentError("Wrong occmat size"))
  end

  out = Vector{Int}(undef, length(phs.sites))
  for (isite, site) in enumerate(phs.sites)
    siteoccbin = compress(PS, occmat[:,isite], BR)
    idx = occupancy_binary_to_state(site, siteoccbin)
    out[isite] = idx
  end
  return out
end


function occupancy_array_to_location_array(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occmat::Matrix{<:Integer}) where {PS, BR, QN}
  n_sites = length(phs.sites)
  n_particles = num_particle_species(PS)
  out = [Int[] for i in 1:n_particles]
  for iptl in 1:n_particles
    for isite in 1:n_sites
      append!(out[iptl], isite for c in 1:occmat[iptl, isite])
    end
  end
  return out
end


function occupancy_binary_to_state_array(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occbin::Unsigned)::Vector{Int} where {PS, BR, QN}
  out = Vector{Int}(undef, length(phs.sites))
  for (isite, site) in enumerate(phs.sites)
    siteoccbin = (occbin & get_bitmask(phs, :, isite)) >> phs.bitoffsets[isite]
    idx = occupancy_binary_to_state(site, siteoccbin)
    out[isite] = idx
  end
  return out
end


function occupancy_binary_to_occupancy_array(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occbin::Unsigned) where {PS, BR, QN}
  n_particles = num_particle_species(PS)
  n_sites = length(phs.sites)
  occmat = Matrix{Int}(undef, (n_particles, n_sites))
  for isite in 1:n_sites
    siteoccbin = (occbin & get_bitmask(phs, :, isite)) >> phs.bitoffsets[isite]
    occmat[:, isite] = extract(PS, siteoccbin)
  end
  return occmat
end


function occupancy_binary_to_location_array(
    phs::ParticleHilbertSpace{PS, BR, QN},
    occbin::Unsigned) where {PS, BR, QN}
  occmat = occupancy_binary_to_occupancy_array(phs, occbin)
  return occupancy_array_to_location_array(phs, occmat)
end


function location_array_to_occupancy_array(
    phs::ParticleHilbertSpace{PS, BR, QN},
    particle_locations::AbstractVector{<:AbstractVector{<:Integer}}) where {PS, BR, QN}

  n_sites = length(phs.sites)
  n_particles = num_particle_species(PS)

  if length(particle_locations) != n_particles
    throw(ArgumentError("length of particle locations does not match the number of particles"))
  end

  sgn = 1
  occmat = zeros(Int, (n_particles, n_sites)) # occupation matrix
  for (iptl, p) in enumerate(particle_locations)
    if exchangesign(particle_species(PS, iptl)) != 1
      sgn *= isparityodd(p) ? -1 : 1
    end
    for isite in p
      occmat[iptl, isite] += 1
    end
  end
  return (occmat, sgn)
end



function location_array_to_occupancy_binary(
        phs::ParticleHilbertSpace{PS, BR, QN},
        particle_locations::AbstractVector{<:AbstractVector{<:Integer}}) where {PS, BR, QN}
  occmat, sgn = location_array_to_occupancy_array(phs, particle_locations)
  return (occupancy_array_to_occupancy_binary(phs, occmat), sgn)
end


"""
    location_array_to_state_array


particles : Vector of (Vector of particle location)

# Example

```
  P S 1 2 3 4 5 = |0, e↓, e↑, m, m⟩
  e↑  0 0 1 0 0
  e↓  0 1 0 0 0 = c†(m,4) c†(m,5) c†(e↑,3) c†(e↓,2) |Ω⟩
  m   0 0 0 1 1
```
"""
function location_array_to_state_array(
        phs::ParticleHilbertSpace{PS, BR, QN},
        particle_locations::AbstractVector{<:AbstractVector{<:Integer}}) where {PS, BR, QN}
  occmat, sgn = location_array_to_occupancy_array(phs, particle_locations)
  return (occupancy_array_to_state_array(phs, occmat), sgn)
end
