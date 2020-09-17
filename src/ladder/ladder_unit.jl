export LadderUnitOperator
export LadderType, CREATION, ANNIHILATION

import LinearAlgebra

@enum LadderType CREATION ANNIHILATION

struct LadderUnitOperator{PS<:ParticleSector, PI, OI}<:AbstractParticleLadderOperator{PS}
    particle_index::PI   # which particle
    orbital::OI          # which orbital
    ladder::LadderType   # creation or annihilation

    function LadderUnitOperator(::Type{PS}, p::P, o::O, l::LadderType) where {PS<:ParticleSector, P, O}
        return new{PS, P, O}(p, o, l)
    end
    function LadderUnitOperator(::PS, p::P, o::O, l::LadderType) where {PS<:ParticleSector, P, O}
        return new{PS, P, O}(p, o, l)
    end
end

function Base.:(==)(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    return (lhs.particle_index == rhs.particle_index) && (lhs.orbital == rhs.orbital) && (lhs.ladder == rhs.ladder)
end


# local normal ordering
function isless_localnormalorder(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    lhs.orbital != rhs.orbital && return isless(lhs.orbital, rhs.orbital)
    lhs.particle_index != rhs.particle_index && return isless(lhs.particle_index, rhs.particle_index)
    return isless(lhs.ladder, rhs.ladder)
end

function isless_normalorder(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    (lhs.ladder != rhs.ladder) && return isless(lhs.ladder, rhs.ladder)
    if lhs.ladder == CREATION # both CREATION
        (lhs.particle_index != rhs.particle_index) && return isless(lhs.particle_index, rhs.particle_index)
        return isless(lhs.orbital, rhs.orbital)
    else # both ANNIHILATION
        (lhs.particle_index != rhs.particle_index) && return isless(rhs.particle_index, lhs.particle_index)
        return isless(rhs.orbital, lhs.orbital)
    end
end


# Particle first
# function Base.isless(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
#     (lhs.ladder != rhs.ladder) && return isless(lhs.ladder, rhs.ladder)
#     if lhs.ladder == CREATION # both CREATION
#         (lhs.particle_index != rhs.particle_index) && return isless(lhs.particle_index, rhs.particle_index)
#         return isless(lhs.orbital, rhs.orbital)
#     else # both ANNIHILATION
#         (lhs.particle_index != rhs.particle_index) && return isless(rhs.particle_index, lhs.particle_index)
#         return isless(rhs.orbital, lhs.orbital)
#     end
# end

# Site first
function Base.isless(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    (lhs.ladder != rhs.ladder) && return isless(lhs.ladder, rhs.ladder)
    if lhs.ladder == CREATION # both CREATION
        (lhs.orbital != rhs.orbital) && return isless(lhs.orbital, rhs.orbital)
        return isless(lhs.particle_index, rhs.particle_index)
    else # both ANNIHILATION
        (lhs.orbital != rhs.orbital) && return isless(rhs.orbital, lhs.orbital)
        return isless(rhs.particle_index, lhs.particle_index)
    end
end

function exchangesign(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    lhs.particle_index != rhs.particle_index && return 1
    isfermion(getspecies(PS, lhs.particle_index)) && return -1
    return 1
end

function maxoccupancy(arg::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    return maxoccupancy(getspecies(PS, arg.particle_index))
end


Base.iszero(arg::LadderUnitOperator) = false

function Base.adjoint(arg::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    new_ladder = arg.ladder == CREATION ? ANNIHILATION : CREATION
    return LadderUnitOperator(PS, arg.particle_index, arg.orbital, new_ladder)
end

LinearAlgebra.ishermitian(arg::LadderUnitOperator) = false




#=

import ExactDiagonalization.apply
import ExactDiagonalization.apply!

function apply!(
    out::SparseState{S1, BR},
    pureop::ParticleProjectorUnitOperator{BR, S2},
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
