export AbstractParticleOperator
export AbstractParticleLadder

abstract type AbstractParticleLadder{PS<:ParticleSector, S<:Number} end
abstract type AbstractParticleOperator{PS<:ParticleSector, S<:Number} <: AbstractOperator{S} end
