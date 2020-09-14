export AbstractParticleOperator
export AbstractParticleLadderOperator
export ParticleCreationOperator
export ParticleAnnihilationOperator

export shortstring


abstract type AbstractParticleOperator{P<:ParticleSector} end
abstract type AbstractParticleLadderOperator{P<:ParticleSector} <: AbstractParticleOperator{P} end

struct ParticleCreationOperator{P<:ParticleSector}<:AbstractParticleLadderOperator{P}
    particleindex::Int
    siteindex::Int

    function ParticleCreationOperator{P}(p::Integer, s::Integer) where {P<:ParticleSector}
        !(1 <= p <= numparticletype(P)) && throw(ArgumentError("particleindex must be within [1, numparticletypes(P)]"))
        return new{P}(p, s)
    end
end

struct ParticleAnnihilationOperator{P<:ParticleSector}<:AbstractParticleLadderOperator{P}
    particleindex::Int
    siteindex::Int

    function ParticleAnnihilationOperator{P}(p::Integer, s::Integer) where {P<:ParticleSector}
        !(1 <= p <= numparticletype(P)) && throw(ArgumentError("particleindex must be within [1, numparticletypes(P)]"))
        return new{P}(p, s)
    end
end

function shortstring(op::ParticleCreationOperator{ParticleSector{P}}) where P
    return "ψ†<sub>$(P.parameters[op.particleindex].parameters[1])</sub>($(op.siteindex))"
end

function shortstring(op::ParticleAnnihilationOperator{ParticleSector{P}}) where P
    return "ψ<sub>$(P.parameters[op.particleindex].parameters[1])</sub>($(op.siteindex))"
end

export apply
function apply(
    phs::ParticleHilbertSpace{P, QN},
    op::ParticleCreationOperator{P},
    bvec::BR
) where {P, QN, BR}

    occ = getparticleoccupation(phs, extract(phs, bvec))
    pi, si = op.particleindex, op.siteindex
    occ[pi, si] += 1

    getstate

end



# @enum LadderType CREATION ANNIHILATION

# struct ParticleLadderOperator{P<:ParticleSector}<:AbstractParticleLadderOperator
#     particleindex::Int
#     siteindex::Int
#     ladder::LadderType

#     function ParticleLadderOperator{P}(p::Integer, s::Integer, l::LadderType) where {P<:ParticleSector}
#         !(1 <= p <= numparticletype(P)) && throw(ArgumentError("particleindex must be within [1, numparticletypes(P)]"))
#         return new{P}(p, s, l)
#     end
# end
