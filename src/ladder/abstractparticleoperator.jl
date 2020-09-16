export AbstractParticleOperator
export AbstractParticleLadderOperator

abstract type AbstractParticleOperator <: AbstractOperator end
abstract type AbstractParticleLadderOperator{PS<:ParticleSector} <: AbstractParticleOperator end
