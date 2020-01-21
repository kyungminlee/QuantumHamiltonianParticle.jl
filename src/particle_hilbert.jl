using ExactDiagonalization
import ExactDiagonalization.get_bitmask
export get_bitmask_particle

export ParticleHilbertSpace

# Add a decoration to the existing Hilbert space
struct ParticleHilbertSpace{PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}
  sites::Vector{ParticleSite{PS, BR, QN}}
  bitwidths::Vector{Int}
  bitoffsets::Vector{Int}

  function ParticleHilbertSpace(sites::AbstractVector{ParticleSite{PS, BR, QN}}) where {PS, BR, QN}
    bitwidths = map(bitwidth, sites)
    bitoffsets = Int[0, cumsum(bitwidths)...]
    if sizeof(BR) * 8 < bitoffsets[end]
      throw(ArgumentError("type $BR too small to represent the hilbert space (need $(bitoffset[end]) bits)"))
    end
    return new{PS, BR, QN}(sites, bitwidths, bitoffsets)
  end
end

bitwidth(hs::ParticleHilbertSpace) = hs.bitoffsets[end]
function bitoffset(phs::ParticleHilbertSpace{PS, BR, QN}, iptl::Integer, isite::Integer) where {PS, BR, QN}
  return phs.bitoffsets[isite] + bitoffset(PS, iptl)
end


function get_bitmask(phs::ParticleHilbertSpace{PS, BR, QN}, iptl::Integer, isite::Integer)::BR where {PS, BR, QN}
  n_sites = length(phs.sites)
  bm = get_bitmask(PS, iptl, BR)
  return bm << phs.bitoffsets[isite]
end

function get_bitmask(phs::ParticleHilbertSpace{PS, BR, QN},
                     iptls::AbstractVector{<:Integer},
                     isites::AbstractVector{<:Integer})::BR where {PS, BR, QN}
  #n_sites = length(phs.sites)
  bm = mapreduce((iptl) -> get_bitmask(PS, iptl, BR), |, zero(BR), iptls)
  return mapreduce((isite) -> bm << phs.bitoffsets[isite], | zero(BR), isites)
end

function get_bitmask(phs::ParticleHilbertSpace{PS, BR, QN}, iptl::Integer, ::Colon)::BR where {PS, BR, QN}
  n_sites = length(phs.sites)
  bm = get_bitmask(PS, iptl, BR)
  return mapreduce( (isite) -> bm << phs.bitoffsets[isite], |, zero(BR), 1:n_sites)
  # out = zero(BR)
  # for isite in 1:n_sites
  #   out |= bm << phs.bitoffsets[isite]
  # end
  # return out
end

function get_bitmask(hs::ParticleHilbertSpace{PS, BR, QN}, ::Colon, isite::Integer)::BR where {PS, BR, QN}
  return make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite], BR)
end


function get_bitmask(hs::ParticleHilbertSpace{PS, BR, QN}, ::Colon, ::Colon)::BR where {PS, BR, QN}
  return make_bitmask(bitwidth(hs), BR)
end


function get_bitmask(hs::ParticleHilbertSpace{PS, BR, QN})::BR where {PS, BR, QN}
  return make_bitmask(bitwidth(hs), BR)
end





function get_bitmask(phs::ParticleHilbertSpace{PS, BR, QN}, iptl::Integer, isites::AbstractVector{<:Integer})::BR where {PS, BR, QN}
  n_sites = length(phs.sites)
  bm = get_bitmask(PS, iptl, BR)
  out = zero(BR)
  for isite in isites
    out |= bm << phs.bitoffsets[isite]
  end
  return out
end




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
