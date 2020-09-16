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

    function LadderSumOperator(op::LadderProductOperator{PS, P, O}, am::S) where {PS, P, O, S<:Number}
        return new{PS, P, O, S}([op=>am,])
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
    return LadderSumOperator([LadderProductOperator(LadderUnitOperator{PS, P, O})=>one(S)])
end

function Base.zero(::Type{LadderSumOperator{PS, P, O, S}}) where {PS, P, O, S}
    return LadderSumOperator(Pair{LadderProductOperator{PS, P, O}, S}[])
end

Base.iszero(arg::LadderSumOperator) = Base.isempty(arg.terms)



function Base.convert(::Type{LadderSumOperator{PS, P, O, S}}, obj::LadderUnitOperator{PS, P, O}) where {PS, P, O, S}
    return LadderSumOperator([LadderProductOperator([obj])=>one(S)])
end

function Base.convert(::Type{LadderSumOperator{PS, P, O, S}}, obj::LadderProductOperator{PS, P, O}) where {PS, P, O, S}
    return LadderSumOperator([obj=>one(S)])
end

# 1. Unary

Base.:(+)(arg::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S} = arg
Base.:(-)(arg::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S} = LadderSumOperator([t => -a for (t, a) in arg.terms])


# 2. Products

# 2a. Simple
Base.:(*)(lhs::LadderUnitOperator, rhs::Number) = LadderSumOperator([LadderProductOperator([lhs]) => rhs])
Base.:(*)(lhs::Number, rhs::LadderUnitOperator) = LadderSumOperator([LadderProductOperator([rhs]) => lhs])
Base.:(/)(lhs::LadderUnitOperator, rhs::Number) = LadderSumOperator([LadderProductOperator([lhs]) => inv(rhs)])
Base.:(\)(lhs::Number, rhs::LadderUnitOperator) = LadderSumOperator([LadderProductOperator([rhs]) => inv(lhs)])
Base.:(//)(lhs::LadderUnitOperator, rhs::Number) = LadderSumOperator([LadderProductOperator([lhs]) => 1//rhs])

Base.:(*)(lhs::LadderProductOperator, rhs::Number) = LadderSumOperator([lhs => rhs])
Base.:(*)(lhs::Number, rhs::LadderProductOperator) = LadderSumOperator([rhs => lhs])
Base.:(/)(lhs::LadderProductOperator, rhs::Number) = LadderSumOperator([lhs => inv(rhs)])
Base.:(\)(lhs::Number, rhs::LadderProductOperator) = LadderSumOperator([rhs => inv(lhs)])
Base.:(//)(lhs::LadderProductOperator, rhs::Number) = LadderSumOperator([lhs => 1//rhs])

Base.:(*)(lhs::LadderSumOperator, rhs::Number) = LadderSumOperator([(t => a * rhs) for (t, a) in lhs.terms])
Base.:(*)(lhs::Number, rhs::LadderSumOperator) = LadderSumOperator([(t => lhs * a) for (t, a) in rhs.terms])
Base.:(/)(lhs::LadderSumOperator, rhs::Number) = LadderSumOperator([(t => a / rhs) for (t, a) in lhs.terms])
Base.:(\)(lhs::Number, rhs::LadderSumOperator) = LadderSumOperator([(t => lhs \ a) for (t, a) in rhs.terms])
Base.:(//)(lhs::LadderSumOperator, rhs::Number) = LadderSumOperator([(t => a // rhs) for (t, a) in lhs.terms])


# 2b. Complex
function Base.:(*)(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return LadderSumOperator([lhs * t => a for (t, a) in rhs.terms])
end

function Base.:(*)(lhs::LadderSumOperator{PS, P, O, S}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O, S}
    return LadderSumOperator([t * rhs => a for (t, a) in lhs.terms])
end

function Base.:(*)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return LadderSumOperator([lhs * t => a for (t, a) in rhs.terms])
end

function Base.:(*)(lhs::LadderSumOperator{PS, P, O, S}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O, S}
    return LadderSumOperator([t * rhs => a for (t, a) in lhs.terms])
end

function Base.:(*)(lhs::LadderSumOperator{PS, P, O, S1}, rhs::LadderSumOperator{PS, P, O, S2}) where {PS, P, O, S1, S2}
    S3 = promote_type(S1, S2)
    terms = Vector{Pair{LadderProductOperator{PS, P, O}, S3}}(undef, length(lhs.terms) * length(rhs.terms))
    i = 0
    for (tl, al) in lhs.terms, (tr, ar) in rhs.terms
        i += 1
        terms[i] = (LadderProductOperator(vcat(tl.factors, tr.factors)) => al * ar)
    end
    return LadderSumOperator(terms)
end


# 3. Addition

function Base.:(+)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return LadderSumOperator([lhs=>one(Int), rhs=>one(Int)])
end

function Base.:(-)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return LadderSumOperator([lhs=>one(Int), rhs=>-one(Int)])
end

Base.:(+)(lhs::LadderSumOperator{PS, P, O, S1}, rhs::LadderSumOperator{PS, P, O, S2}) where {PS, P, O, S1, S2} = LadderSumOperator(vcat(lhs.terms, rhs.terms))
# Base.:(-)(lhs::LadderSumOperator{PS, P, O, S1}, rhs::LadderSumOperator{PS, P, O, S2}) where {PS, P, O, S1, S2} = (lhs) + (-rhs)

function Base.:(+)(lhs::LadderSumOperator{PS, P, O, S}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O, S}
    return LadderSumOperator([lhs.terms..., rhs=>one(S)])
end

function Base.:(+)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return LadderSumOperator([lhs=>one(S), rhs.terms...])
end

function isequiv(lhs::LadderSumOperator{PS, P, O, S1}, rhs::LadderSumOperator{PS, P, O, S2}) where {PS, P, O, S1, S2}
    return isempty(simplify(lhs - rhs))
end


function Base.adjoint(arg::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return LadderSumOperator([(adjoint(t), conj(a)) for (t, a) in arg.terms])
end

function LinearAlgebra.ishermitian(arg::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return isequiv(arg, adjoint(arg))
end
