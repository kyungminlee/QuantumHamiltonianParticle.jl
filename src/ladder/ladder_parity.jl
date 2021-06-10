export get_fermion_parity


"""
    get_fermion_parity(phs, op, bvec)

Get the fermion parity (0 or 1) for when op is applied to bvec.
"""
function get_fermion_parity(
    hs::ParticleHilbertSpace,
    op::ParticleLadderUnit{PS, <:Integer, <:Integer},
    bvec::Unsigned,
) where {PS}
    bm_mask = get_parity_bitmask(hs, op.particle_index, op.orbital)
    bm_parity = bm_mask & bvec
    return count_ones(bm_parity) % 2
end



"""
    get_fermion_parity(phs, op, bvec)

Get the fermion parity (0 or 1) for when op is applied to bvec.
"""
function get_fermion_parity(
    op::ParticleLadderUnit{PS, <:Integer, <:Integer},
    bvec::Unsigned,
) where {PS}
    bm_mask = get_parity_bitmask(PS, op.particle_index, op.orbital)
    bm_parity = bm_mask & bvec
    return count_ones(bm_parity) % 2
end
