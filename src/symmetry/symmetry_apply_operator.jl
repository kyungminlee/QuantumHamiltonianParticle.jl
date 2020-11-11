using LatticeTools

# import ExactDiagonalization.AbstractSymmetryOperation
import ExactDiagonalization.symmetry_apply

# === With Hilbert Space ===
function symmetry_apply(
    ::ParticleHilbertSpace{PS, BR, QN},
    permutation::SitePermutation,
    op::AbstractParticleLadder{PS},
) where {PS, BR, QN}
    return symmetry_apply(permutation, op)
end

# === Ladders ===
function symmetry_apply(permutation::SitePermutation, op::ParticleLadderUnit{PS, <:Integer, <:Integer}) where PS
    return ParticleLadderUnit(PS, op.particle_index, permutation(op.orbital), op.ladder)
end

# No need to worry about signs here. The operators are not ordered.
function symmetry_apply(permutation::SitePermutation, op::ParticleLadderProduct{PI, OI}) where {PI, OI}
    return ParticleLadderProduct([symmetry_apply(permutation, f) for f in op.factors])
end

function symmetry_apply(permutation::SitePermutation, op::ParticleLadderSum{PI, OI}) where {PI, OI}
    return ParticleLadderSum([symmetry_apply(permutation, t) => a for (t, a) in op.terms])
end

# === Ladder Embedding ===
function symmetry_apply(permutation::SitePermutation, op::ParticleLadderOperatorEmbedding)
    return ParticleLadderOperatorEmbedding(op.hilbert_space, symmetry_apply(permutation, op.ladder))
end


import ExactDiagonalization.isinvariant

function isinvariant(permutation::SitePermutation, op::AbstractParticleOperator)
    return iszero(simplify(symmetry_apply(permutation, op) - op))
end

function isinvariant(permutation::SitePermutation, op::AbstractParticleLadder)
    return iszero(simplify(symmetry_apply(permutation, op) - op))
end
