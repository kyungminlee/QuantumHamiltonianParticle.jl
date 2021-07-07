export ParticleLadderUnit
export LadderType, CREATION, ANNIHILATION

import LinearAlgebra

@enum LadderType CREATION ANNIHILATION

struct ParticleLadderUnit{PS<:ParticleSector, PI, OI}<:AbstractParticleLadder{PS, Int}
    particle_index::PI   # which particle
    orbital::OI          # which orbital
    ladder::LadderType   # creation or annihilation

    function ParticleLadderUnit(::Type{PS}, p::P, o::O, l::LadderType) where {PS<:ParticleSector, P, O}
        return new{PS, P, O}(p, o, l)
    end
    function ParticleLadderUnit(::PS, p::P, o::O, l::LadderType) where {PS<:ParticleSector, P, O}
        return new{PS, P, O}(p, o, l)
    end
end

function Base.:(==)(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
    return (lhs.particle_index == rhs.particle_index) && (lhs.orbital == rhs.orbital) && (lhs.ladder == rhs.ladder)
end


# # local normal ordering
# function Base.isless(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
#     lhs.orbital != rhs.orbital && return isless(lhs.orbital, rhs.orbital)
#     lhs.particle_index != rhs.particle_index && return isless(lhs.particle_index, rhs.particle_index)
#     return isless(lhs.ladder, rhs.ladder)
# end


# normal ordering
function Base.isless(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
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
# function Base.isless(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
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
# function Base.isless(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
#     (lhs.ladder != rhs.ladder) && return isless(lhs.ladder, rhs.ladder)
#     if lhs.ladder == CREATION # both CREATION
#         (lhs.orbital != rhs.orbital) && return isless(lhs.orbital, rhs.orbital)
#         return isless(lhs.particle_index, rhs.particle_index)
#     else # both ANNIHILATION
#         (lhs.orbital != rhs.orbital) && return isless(rhs.orbital, lhs.orbital)
#         return isless(rhs.particle_index, lhs.particle_index)
#     end
# end

function exchangesign(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
    lhs.particle_index != rhs.particle_index && return 1
    isfermion(getspecies(PS, lhs.particle_index)) && return -1
    return 1
end

function maxoccupancy(arg::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
    return maxoccupancy(getspecies(PS, arg.particle_index))
end


Base.iszero(arg::ParticleLadderUnit) = false

function Base.adjoint(arg::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
    new_ladder = arg.ladder == CREATION ? ANNIHILATION : CREATION
    return ParticleLadderUnit(PS, arg.particle_index, arg.orbital, new_ladder)
end

LinearAlgebra.ishermitian(arg::ParticleLadderUnit) = false
