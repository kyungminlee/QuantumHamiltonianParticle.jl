using LatticeTools
import QuantumHamiltonian.symmetry_apply

function symmetry_apply(
    phs::ParticleHilbertSpace,
    permutation::SitePermutation,
    occbin::BR,
    phase::S=one(Int)
) where {BR<:Unsigned, S<:Number}
    lv = occbin2locvec(phs, occbin)
    lv2 = [permutation.(x) for x in lv] # per particle
    occbin2, sgn = locvec2occbin(phs, lv2)
    return (occbin2, sgn*phase)
end

function symmetry_apply(
    phs::ParticleHilbertSpace,
    perm::LocalGeneralizedPermutation{A},
    occbin::BR,
    amplitude::S=one(Int)
) where {A, BR<:Unsigned, S<:Number}
    @boundscheck let
        if length(phs.sites) != length(perm.operations)
            throw(ArgumentError("number of sites should match"))
        end
        for (isite, (site, op)) in enumerate(zip(phs.sites, perm.operations))
            if dimension(site) != length(op.map)
                throw(ArgumentError("dimension mismatch at site $isite: $(dimension(site)) != $(length(op.map))"))
            end
        end
    end
    sv = occbin2statevec(phs, occbin)
    sv2 = [
        let # this newline is relevant
            y, amplitude = op(x, amplitude)
            y
        end
        for (x, op) in zip(sv, perm.operations)        
    ]
    ob2 = statevec2occbin(phs, sv2)
    return (ob2, amplitude)
end
