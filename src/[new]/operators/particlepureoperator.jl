export ParticlePureOperator


# Convention:
# wave function
#   c†(1) c†(2) c†(3) |0>
#   Note: Wigner-Jordan:
#
#   |0 0| |1  0| = |0 0|
#   |1 0| |0 -1|   |1 0|
#
#   c†(3)      111111-zz
#   c†(7)      11-zzzzzz
# ---------- = ---------
# c†(3)c†(7)   11-zzz-11
#
#   |1  0| |0 0| = | 0 0|
#   |0 -1| |1 0|   |-1 0|
#
#   c†(7)        (11-zzzzzz)
#   c†(3)        (111111-zz)
# ----------  =   ---------
# c†(7)c†(3)    -(11-zzz-11) = -c†(3)c†(7)
#
#   c(3)    11111
#
#
#



struct ParticlePureOperator{P<:ParticleSector, BR<:Unsigned, S<:Number}<:AbstractParticleOperator{P}
    bitmask::BR
    bitrow::BR
    bitcol::BR

    #pbitmask::BR  # basically a duplicate.
    plocrow::Vector{Vector{Int}}
    ploccol::Vector{Vector{Int}}

    amplitude::S
end


function particlecount(hs::ParticleHilbertSpace{P}, pbitvec::BR, loc::Integer) where {P}


end

function popparticle!(loclist::AbstractVector{<:Integer}, loc::Integer)
    # assume loclist is sorted

end





function apply(hs::ParticleHilbertSpace{QN, P}, op::ParticlePureOperator{P}, bvec::Unsigned) {QN, P}
    locvec = getparticlelocationlist(hs, bvec)
    for (lc, lv) in zip(op.loccol, locvec)
        # remove every element of lc from lv
        for l in lc
            delete!

    end
end
