using LatticeTools
import QuantumHamiltonian.symmetry_apply

function symmetry_apply_naive(
    phs::ParticleHilbertSpace,
    permutation::SitePermutation,
    occbin::BR,
    phase::S=one(Int)
) where {BR<:Unsigned, S<:Number}
    n_particles = speciescount(phs)
    n_sites = numsites(phs)
    lv = [Int[] for _ in 1:n_particles]
    for iptl in 1:n_particles
        sizehint!(lv[iptl], n_sites)
    end
    occbin2locvec!(lv, phs, occbin)
    # lv = occbin2locvec(phs, occbin)
    for iptl in 1:n_particles
        for j in 1:length(lv[iptl])
            lv[iptl][j] = permutation(lv[iptl][j])
        end
    end
    # lv2 = [permutation.(x) for x in lv] # per particle
    occbin2, sgn = locvec2occbin(phs, lv)
    return (occbin2, sgn*phase)
end

function symmetry_apply(
    phs::ParticleHilbertSpace{PS, BR, QN},
    permutation::SitePermutation,
    occbin::Unsigned,
    phase::S=one(Int)
) where {PS, BR<:Unsigned, QN, S<:Number}
    n_particles::Int = speciescount(phs)
    n_sites::Int = numsites(phs)
    occbin2 = zero(BR)
    pshift = 0
    for iptl in 1:n_particles
        # pbw::Int = bitwidth(PS, iptl) # Keep this line for debugging purpose
        # pbm::BR = make_bitmask(pbw, BR) # Keep this line for debugging purpose
        pbw::Int = bitwidth(phs, iptl, 1)
        pbm::BR = get_bitmask(phs, iptl, 1) >> bitoffset(phs, iptl, 1)
        p_occbin2 = zero(BR)
        if exchangesign(PS, iptl) == 1
            for isite in 1:n_sites
                jsite = permutation(isite)
                occ = (occbin >> bitoffset(phs, isite)) & pbm
                p_occbin2 |= occ << bitoffset(phs, jsite)
            end
        elseif exchangesign(PS, iptl) == -1
            for isite in 1:n_sites
                jsite = permutation(isite)
                occ = (occbin >> bitoffset(phs, isite)) & pbm
                p_occbin2 |= occ << bitoffset(phs, jsite)
                if occ % 2 != 0
                    count = 0
                    for ksite in jsite+1:n_sites
                        count += Int((p_occbin2 >> bitoffset(phs, ksite)) & pbm)
                    end
                    if count % 2 != 0
                        phase = -phase
                    end
                end
            end
        else
            throw(ArgumentError("Unsupported exchange sign $(exchangesign(getspecies(PS, iptl))) for $(getspecies(PS, iptl))")) # COV_EXCL_LINE
        end
        occbin >>= pbw
        occbin2 |= p_occbin2 << pshift
        pshift += pbw
    end
    return (occbin2, phase)
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
    # the semicolon is relevant
    sv2 = [let ; y, amplitude = op(x, amplitude); y end
           for (x, op) in zip(sv, perm.operations)]
    ob2 = statevec2occbin(phs, sv2)
    return (ob2, amplitude)
end
