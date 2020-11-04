export get_fermion_parity


function get_fermion_parity(
    hs::ParticleHilbertSpace,
    op::ParticleLadderUnit{PS, <:Integer, <:Integer},
    bvec::Unsigned,
) where {PS}
    bm_mask = get_parity_bitmask(hs, op.particle_index, op.orbital)
    bm_parity = bm_mask & bvec
    return count_ones(bm_parity) % 2
end
