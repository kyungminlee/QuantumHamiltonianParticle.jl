export OrbitalPermutation

using GroupTools

export OrbitalPermutation
export orbitalpermutation

struct OrbitalPermutation{A} <: GroupTools.AbstractSymmetryOperation
    operation::GeneralizedPermutation{A}
    OrbitalPermutation(operations::GeneralizedPermutation{A}) where {A} = new{A}(operations)
    OrbitalPermutation(map::AbstractVector{<:Integer}, phase::AbstractVector{<:Phase}) = OrbitalPermutation(GeneralizedPermutation(map, phase))
end

function Base.:(*)(lhs::OrbitalPermutation, rhs::OrbitalPermutation)
    return OrbitalPermutation(lhs.operation * rhs.operation)
end

function LatticeTools.isidentity(arg::OrbitalPermutation)
    return all(LatticeTools.isidentity, arg.operation)
end

Base.inv(arg::OrbitalPermutation) = OrbitalPermutation(inv(arg.operation))

Base.hash(arg::L, h::UInt) where {L<:OrbitalPermutation} = hash(L, hash(arg.operation, h))

function Base.:(==)(lhs::OrbitalPermutation, rhs::OrbitalPermutation)
    return all(lhs.operation == rhs.operation)
end

(permutation::SitePermutation)(op::OrbitalPermutation) = op

function symmetry_apply(
    hs::ParticleHilbertSpace,
    perm::OrbitalPermutation{A},
    bitrep::BR,
    amplitude::S=one(Int)
) where {A, BR<:Unsigned, S<:Number}
    occmat = occbin2occmat(hs, bitrep)
    (nptl, nsite) = size(occmat)
    newoccmat = zeros(Int, (nptl, nsite))

    for isite in 1:nsite, iptl in 1:nptl
        p = occmat[iptl, isite]
        jptl = perm.operation.map[iptl]
        am = perm.operation.phase[iptl]
        newoccmat[jptl, isite] = p
        amplitude = (am^p)(amplitude)
    end
    newbinrep = occmat2occbin(hs, newoccmat)
    return (newbinrep, amplitude)
end


"""
# Example

```
orbitalpermutation(PS, :px => (:py, 1), :py => (:px, -1))
```
"""
function orbitalpermutation(
    ::Type{PS},
    mapping::Pair{Symbol, <:Tuple{Symbol, <:Phase}}...
) where {PS<:ParticleSector}
    np = numspecies(PS)
    mapping_index = collect(1:np)
    mapping_phase = ones(Phase{Rational{Int}}, np)
    for (s, (t, a)) in mapping
        si = findspeciesindex(PS, s)
        ti = findspeciesindex(PS, t)
        si > 0 || continue
        ti > 0 || continue
        mapping_index[si] = ti
        mapping_phase[si] = a
    end
    return OrbitalPermutation(mapping_index, mapping_phase)
end

orbitalpermutation(::PS, args...) where {PS<:ParticleSector} = orbitalpermutation(PS, args...)