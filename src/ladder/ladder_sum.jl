export ParticleLadderSum

import LinearAlgebra
import ExactDiagonalization.isequiv

struct ParticleLadderSum{PS<:ParticleSector, P, O, S<:Number}<:AbstractParticleLadder{PS, S}
    terms::Vector{Pair{ParticleLadderProduct{PS, P, O}, S}}

    function ParticleLadderSum{PS, P, O, S}() where {PS, P, O, S}
        return new{PS, P, O, S}([])
    end

    function ParticleLadderSum(op::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
        return new{PS, P, O, Int}([ParticleLadderProduct([op])=>1,])
    end

    function ParticleLadderSum(op::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
        return new{PS, P, O, Int}([op=>1,])
    end

    function ParticleLadderSum(ops::Vararg{Pair{ParticleLadderProduct{PS, P, O}, S}}) where {PS, P, O, S<:Number}
        return new{PS, P, O, S}([op=>am for (op, am) in ops])
    end

    function ParticleLadderSum(terms::AbstractVector{Pair{ParticleLadderProduct{PS, P, O}, S}}) where {PS, P, O, S<:Number}
        return new{PS, P, O, S}(terms)
    end
end

Base.:(==)(lhs::OP, rhs::OP) where {OP<:ParticleLadderSum} = lhs.terms == rhs.terms

function Base.one(::Type{ParticleLadderSum{PS, P, O, S}}) where {PS, P, O, S}
    return ParticleLadderSum([one(ParticleLadderProduct{PS, P, O}) => one(S)])
end

function Base.zero(::Type{ParticleLadderSum{PS, P, O, S}}) where {PS, P, O, S}
    return ParticleLadderSum(Pair{ParticleLadderProduct{PS, P, O}, S}[])
end

Base.one(::OP) where {OP<:ParticleLadderSum} = Base.one(OP)
Base.zero(::OP) where {OP<:ParticleLadderSum} = Base.zero(OP)

Base.iszero(arg::ParticleLadderSum) = Base.isempty(arg.terms)


function isequiv(lhs::ParticleLadderSum{PS, P, O, S1}, rhs::ParticleLadderSum{PS, P, O, S2}) where {PS, P, O, S1, S2}
    return iszero(simplify(lhs - rhs))
end


function Base.adjoint(arg::ParticleLadderSum{PS, P, O, S}) where {PS, P, O, S}
    return ParticleLadderSum([(adjoint(t)=>conj(a)) for (t, a) in arg.terms])
end

function LinearAlgebra.ishermitian(arg::ParticleLadderSum{PS, P, O, S}) where {PS, P, O, S}
    return isequiv(arg, adjoint(arg))
end
