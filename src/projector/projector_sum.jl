export ParticleProjectorSumOperator

struct ParticleProjectorSumOperator{BR<:Unsigned, S<:Number}<:AbstractParticleProjectorOperator{BR, S}
    terms::Vector{ParticleProjectorUnitOperator{BR, S}}
    function ParticleProjectorSumOperator(terms::AbstractVector{ParticleProjectorUnitOperator{BR, S}}) where {BR, S}
        return new{BR, S}(terms)
    end
end


function Base.zero(::Type{ParticleProjectorSumOperator{BR, S}}) where {BR, S}
    return ParticleProjectorSumOperator(ParticleProjectorUnitOperator{BR, S}[])
end

function Base.one(::Type{ParticleProjectorSumOperator{BR, S}}) where {BR, S}
    return ParticleProjectorSumOperator([one(ParticleProjectorUnitOperator{BR, S})])
end

Base.iszero(arg::ParticleProjectorSumOperator) = isempty(arg.terms)


# Equality

function Base.:(==)(x::ParticleProjectorSumOperator, y::ParticleProjectorSumOperator)
    return x.terms == y.terms
end

# function Base.isapprox(
#     x::ParticleProjectorSumOperator, y::ParticleProjectorSumOperator;
#     atol::Real=0, rtol::Real=Base.rtoldefault(x.amplitude,y.amplitude,atol), nans::Bool=false,
# )
#     return isapprox(x.terms, y.terms; atol=atol, rtol=rtol, nans=nans)
# end


# Unary

function Base.real(x::ParticleProjectorSumOperator)
    return ParticleProjectorSumOperator(real.(x.terms))
end

function Base.imag(x::ParticleProjectorSumOperator)
    return ParticleProjectorSumOperator(imag.(x.terms))
end

function Base.adjoint(x::ParticleProjectorSumOperator)
    return ParticleProjectorSumOperator(adjoint.(x.terms))
end

function Base.conj(x::ParticleProjectorSumOperator)
    return ParticleProjectorSumOperator(conj.(x.terms))
end

function Base.transpose(x::ParticleProjectorSumOperator)
    return ParticleProjectorSumOperator(transpose.(x.terms))
end

Base.:(+)(arg::ParticleProjectorSumOperator) = arg
Base.:(-)(arg::ParticleProjectorSumOperator) = ParticleProjectorSumOperator(-arg.terms)


# Binary (scale)

Base.:(*)(x::ParticleProjectorSumOperator, y::Number)  = ParticleProjectorSumOperator(x.terms .* y)
Base.:(/)(x::ParticleProjectorSumOperator, y::Number)  = ParticleProjectorSumOperator(x.terms ./ y)
Base.:(//)(x::ParticleProjectorSumOperator, y::Number) = ParticleProjectorSumOperator(x.terms .// y)
Base.:(*)(y::Number, x::ParticleProjectorSumOperator)  = ParticleProjectorSumOperator(y .* x.terms)
Base.:(\)(y::Number, x::ParticleProjectorSumOperator)  = ParticleProjectorSumOperator(y .\ x.terms)

function Base.:(+)(x::ParticleProjectorSumOperator{BR, S1}, y::S2) where {BR, S1, S2}
    S = promote_type(S1, S2)
    return ParticleProjectorSumOperator(
        ParticleProjectorUnitOperator{BR, S}[x.terms..., ParticleProjectorUnitOperator{BR, S}(y)]
    )
end

function Base.:(+)(y::S2, x::ParticleProjectorSumOperator{BR, S1}) where {BR, S1, S2}
    S = promote_type(S1, S2)
    return ParticleProjectorSumOperator(
        ParticleProjectorUnitOperator{BR, S}[ParticleProjectorUnitOperator{BR, S}(y), x.terms...]
    )
end

Base.:(+)(lhs::ParticleProjectorUnitOperator, rhs::ParticleProjectorUnitOperator) = ParticleProjectorSumOperator([lhs, rhs])
Base.:(+)(lhs::ParticleProjectorUnitOperator, rhs::ParticleProjectorSumOperator) = ParticleProjectorSumOperator([lhs, rhs.terms...])
Base.:(+)(lhs::ParticleProjectorSumOperator, rhs::ParticleProjectorUnitOperator) = ParticleProjectorSumOperator([lhs.terms..., rhs])
Base.:(+)(lhs::ParticleProjectorSumOperator, rhs::ParticleProjectorSumOperator) = ParticleProjectorSumOperator(vcat(lhs.terms, rhs.terms))

function Base.:(*)(x::ParticleProjectorSumOperator{BR, S1}, y::ParticleProjectorSumOperator{BR, S2}) where {BR, S1, S2}
    S3 = promote_type(S1, S2)
    terms = ParticleProjectorUnitOperator{BR, S3}[]
    for t1 in x.terms, t2 in y.terms
        t3 = t1 * t2
        !iszero(t3) && push!(terms, t3)
    end
    return ParticleProjectorSumOperator(terms)
end
