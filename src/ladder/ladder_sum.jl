export LadderSumOperator

import LinearAlgebra
import ExactDiagonalization.isequiv

struct LadderSumOperator{PS, P, O, S<:Number}<:AbstractParticleLadderOperator{PS}
    terms::Vector{Pair{LadderProductOperator{PS, P, O}, S}}

    function LadderSumOperator(op::LadderUnitOperator{PS, P, O}) where {PS, P, O}
        return new{PS, P, O, Int}([LadderProductOperator([op])=>1,])
    end

    function LadderSumOperator(op::LadderProductOperator{PS, P, O}) where {PS, P, O}
        return new{PS, P, O, Int}([op=>1,])
    end

    function LadderSumOperator(ops::Vararg{Pair{LadderProductOperator{PS, P, O}, S}}) where {PS, P, O, S<:Number}
        return new{PS, P, O, S}([op=>am for (op, am) in ops])
    end

    function LadderSumOperator(terms::AbstractVector{Pair{LadderProductOperator{PS, P, O}, S}}) where {PS, P, O, S<:Number}
        return new{PS, P, O, S}(terms)
    end
end


function Base.:(==)(
    lhs::LadderSumOperator{PS, P, O, S},
    rhs::LadderSumOperator{PS, P, O, S}
) where {PS, P, O, S}
    return lhs.terms == rhs.terms
end

function Base.one(::Type{LadderSumOperator{PS, P, O, S}}) where {PS, P, O, S}
    return LadderSumOperator([one(LadderProductOperator{PS, P, O}) => one(S)])
end

function Base.zero(::Type{LadderSumOperator{PS, P, O, S}}) where {PS, P, O, S}
    return LadderSumOperator(Pair{LadderProductOperator{PS, P, O}, S}[])
end

Base.one(::OP) where {OP<:LadderSumOperator} = Base.one(OP)
Base.zero(::OP) where {OP<:LadderSumOperator} = Base.zero(OP)

Base.iszero(arg::LadderSumOperator) = Base.isempty(arg.terms)

function Base.convert(::Type{LadderSumOperator{PS, P, O, S}}, obj::LadderUnitOperator{PS, P, O}) where {PS, P, O, S}
    return LadderSumOperator([LadderProductOperator([obj])=>one(S)])
end

function Base.convert(::Type{LadderSumOperator{PS, P, O, S}}, obj::LadderProductOperator{PS, P, O}) where {PS, P, O, S}
    return LadderSumOperator([obj=>one(S)])
end

function Base.convert(::Type{LadderSumOperator{PS, P, O, S}}, obj::LadderSumOperator{PS, P, O, S2}) where {PS, P, O, S}
    return LadderSumOperator([t=>convert(S, a) for (t, a) in obj.terms])
end


function isequiv(lhs::LadderSumOperator{PS, P, O, S1}, rhs::LadderSumOperator{PS, P, O, S2}) where {PS, P, O, S1, S2}
    return iszero(simplify(lhs - rhs))
end


function Base.adjoint(arg::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return LadderSumOperator([(adjoint(t)=>conj(a)) for (t, a) in arg.terms])
end

function LinearAlgebra.ishermitian(arg::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return isequiv(arg, adjoint(arg))
end
