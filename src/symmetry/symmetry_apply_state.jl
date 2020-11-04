using LatticeTools
import ExactDiagonalization.symmetry_apply

function symmetry_apply(
    phs::ParticleHilbertSpace,
    permutation::SitePermutation,
    occbin::BR
) where {BR<:Unsigned}
    lv = occbin2locvec(phs, occbin)
    lv2 = [permutation.(x) for x in lv] # per particle
    occbin2, sgn = locvec2occbin(phs, lv2)
    return (occbin2, sgn)
end
