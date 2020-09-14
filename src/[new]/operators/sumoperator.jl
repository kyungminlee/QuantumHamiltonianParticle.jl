
struct ParticleSumOperator{P<:ParticleSector, Scalar}
    terms::Vector{AbstractParticleOperator{P}}
end
