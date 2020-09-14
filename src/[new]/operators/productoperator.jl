struct ParticleProductOperator{P<:ParticleSector, Scalar}
    factors::Vector{AbstractParticleLadderOperator{P}}
    amplitude::Scalar
end

function Base.:(*)(lhs::AbstractParticleLadderOperator{P}, rhs::AbstractParticleLadderOperator{P}) where {P}
    return ParticleProductOperator{P, Bool}([lhs, rhs], true)
end

function Base.:(*)(lhs::AbstractParticleLadderOperator{P}, rhs::ParticleProductOperator{P, S}) where {P, S}
    return ParticleProductOperator{P, S}([lhs, rhs.factors...], rhs.amplitude)
end

function Base.:(*)(lhs::ParticleProductOperator{P, S}, rhs::AbstractParticleLadderOperator{P}) where {P, S}
    return ParticleProductOperator{P, S}([lhs.factors..., rhs], lhs.amplitude)
end

function Base.:(*)(lhs::ParticleProductOperator{P, S1}, rhs::ParticleProductOperator{P, S2}) where {P, S1, S2}
    S3 = promote_type(S1, S2)
    amplitude = lhs.amplitude * rhs.amplitude
    return ParticleProductOperator{P, S3}([lhs.factors..., rhs.factors...], amplitude)
end


function Base.convert(::Type{ParticleProductOperator{P, S}}, arg::AbstractParticleLadderOperator{P}) where {P, S}
    return ParticleProductOperator{P, S}([arg], one(Scalar))
end




function shortstring(op::ParticleProductOperator)
    return (isone(op.amplitude) ? "" : "$(op.amplitude)")*join(shortstring.(op.factors), "⋅")
end
