export get_fermion_parity

function get_parity_mask(hs::ParticleHilbertSpace, iptl::Integer, isite::Integer)
    bm_species = get_bitmask(hs, iptl, :)
    bm_site = get_bitmask(hs, iptl, isite)
    bm_mask = ( bm_site - 1 ) & bm_species  # σᶻ in jordan wigner string
    return bm_mask
end


function get_fermion_parity(
    hs::ParticleHilbertSpace,
    op::LadderUnitOperator{PS, <:Integer, <:Integer},
    bvec::Unsigned,
) where {PS}
    bm_mask = get_parity_mask(hs, op.particle_index, op.orbital)
    bm_parity = bm_mask & bvec
    return count_ones(bm_parity) % 2
end
