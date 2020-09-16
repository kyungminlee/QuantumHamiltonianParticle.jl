export AbstractParticleOperator
export AbstractParticleLadderOperator

abstract type AbstractParticleOperator{PS<:ParticleSector} <: AbstractOperator{Int} end
abstract type AbstractParticleLadderOperator{PS<:ParticleSector} <: AbstractParticleOperator{PS} end
