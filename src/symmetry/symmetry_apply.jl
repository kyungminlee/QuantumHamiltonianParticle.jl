using LatticeTools
import ExactDiagonalization.symmetry_apply

function symmetry_apply(
    # hs::ParticleHilbertSpace,
    permutation::SitePermutation,
    op::LadderUnitOperator{<:Integer, <:Integer},
)
    return LadderUnitOperator(op.particle_index, permutation(op.orbital), op.ladder)
end


function symmetry_apply(
    # hs::ParticleHilbertSpace,
    permutation::SitePermutation,
    op::LadderProductOperator{PI, OI},
) where {PI, OI}
    return LadderProductOperator(symmetry_apply(hs, permutation, f) for f in op.factors)
end

function symmetry_apply(
    # hs::ParticleHilbertSpace,
    permutation::SitePermutation,
    op::LadderSumOperator{PI, OI},
    
)