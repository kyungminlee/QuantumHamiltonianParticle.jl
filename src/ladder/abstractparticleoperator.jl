export AbstractParticleOperator
export AbstractParticleLadder

abstract type AbstractParticleOperator{PS<:ParticleSector, S<:Number} <: AbstractOperator{S} end
abstract type AbstractParticleLadder{PS<:ParticleSector, S<:Number} <: AbstractParticleOperator{PS, S} end
