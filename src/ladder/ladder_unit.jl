export LadderUnitOperator

import LinearAlgebra


export LadderType, CREATION, ANNIHILATION
@enum LadderType CREATION ANNIHILATION

struct LadderUnitOperator{PI, OI}<:AbstractParticleOperator
    particle_index::PI   # which particle
    orbital::OI          # which orbital
    ladder::LadderType   # creation or annihilation

    function LadderUnitOperator(p::P, o::O, l::LadderType) where {P, O}
        return new{P, O}(p, o, l)
    end
end


function maxoccupancy(arg::LadderUnitOperator{PI, OI}) where {PI, OI}
    return maxoccupancy(particle_species(arg.particle_index))
end


function Base.:(==)(lhs::LadderUnitOperator{P, O}, rhs::LadderUnitOperator{P, O}) where {P, O}
    return lhs.particle_index == rhs.particle_index && lhs.orbital == rhs.orbital && lhs.ladder == rhs.ladder
end


function Base.isless(lhs::LadderUnitOperator{P, O}, rhs::LadderUnitOperator{P, O}) where {P, O}
    lhs.ladder != rhs.ladder && return isless(lhs.ladder, rhs.ladder)

    if lhs.ladder == CREATION # both CREATION
        lhs.particle_index != rhs.particle_index && return isless(lhs.particle_index, rhs.particle_index)
        return isless(lhs.orbital, rhs.orbital)
    else # both ANNIHILATION
        lhs.particle_index != rhs.particle_index && return isless(rhs.particle_index, lhs.particle_index)
        return isless(rhs.orbital, lhs.orbital)
    end
end


function exchangesign(lhs::LadderUnitOperator{P, O}, rhs::LadderUnitOperator{P, O}) where {P, O}
    lhs.particle_index != rhs.particle_index && return 1
    isfermion(particle_species(lhs.particle_index)) && return -1
    return 1
end

Base.iszero(arg::LadderUnitOperator) = false




function Base.adjoint(arg::LadderUnitOperator{P, O}) where {P, O}
    new_ladder = arg.ladder == CREATION ? ANNIHILATION : CREATION
    return LadderUnitOperator(arg.particle_index, arg.orbital, new_ladder)
end

LinearAlgebra.ishermitian(arg::LadderUnitOperator) = false

#=

import ExactDiagonalization.apply
import ExactDiagonalization.apply!

function apply!(
    out::SparseState{S1, BR},
    pureop::ParticleProjectionUnitOperator{BR, S2},
    psi::SparseState{S3, BR},
) where {S1, S2, S3, BR}
    for (b, v) in psi.components
        if (b & pureop.bitmask) == pureop.bitcol
            b2 = (b & ~pureop.bitmask) | pureop.bitrow
            a = pureop.amplitude * v
            isparityeven = mod(count_ones(pureop.parity_bitmask & bvec), 2) == 0
            out[b2] += isparityeven ? a : -a
        end
    end
    return out
end
=#
