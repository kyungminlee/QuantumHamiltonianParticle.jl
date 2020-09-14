export LadderSumOperator

import LinearAlgebra

struct LadderSumOperator{P, O, S<:Number}<:AbstractParticleOperator
    terms::Vector{Pair{LadderProductOperator{P, O}, S}}

    function LadderSumOperator(op::LadderUnitOperator{P, O}) where {P, O}
        return new{P, O, Int}([LadderProductOperator([op])=>1,])
    end

    function LadderSumOperator(op::LadderProductOperator{P, O}) where {P, O}
        return new{P, O, Int}([op=>1,])
    end

    function LadderSumOperator(op::LadderProductOperator{P, O}, am::S) where {P, O, S<:Number}
        return new{P, O, S}([op=>am,])
    end

    function LadderSumOperator(terms::AbstractVector{Pair{LadderProductOperator{P, O}, S}}) where {P, O, S<:Number}
        return new{P, O, S}(terms)
    end
end


function Base.:(==)(lhs::LadderSumOperator{P, O, S}, rhs::LadderSumOperator{P, O, S}) where {P, O, S}
    return lhs.terms == rhs.terms
end

function Base.one(::Type{LadderSumOperator{P, O, S}}) where {P, O, S}
    return LadderSumOperator([LadderProductOperator(LadderUnitOperator{P, O})=>one(S)])
end

function Base.zero(::Type{LadderSumOperator{P, O, S}}) where {P, O, S}
    return LadderSumOperator(Pair{LadderProductOperator{P, O}, S}[])
end

Base.iszero(arg::LadderSumOperator) = Base.isempty(arg.terms)



function Base.convert(::Type{LadderSumOperator{P, O, S}}, obj::LadderUnitOperator{P, O}) where {P, O, S}
    return LadderSumOperator([LadderProductOperator([obj])=>one(S)])
end

function Base.convert(::Type{LadderSumOperator{P, O, S}}, obj::LadderProductOperator{P, O}) where {P, O, S}
    return LadderSumOperator([obj=>one(S)])
end

# 1. Unary

Base.:(+)(arg::LadderSumOperator{P, O, S}) where {P, O, S} = arg
Base.:(-)(arg::LadderSumOperator{P, O, S}) where {P, O, S} = LadderSumOperator([(t, -a) for (t, a) in arg.terms])


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
function Base.:(*)(lhs::LadderUnitOperator{P, O}, rhs::LadderSumOperator{P, O, S}) where {P, O, S}
    return LadderSumOperator([lhs * t => a for (t, a) in rhs.terms])
end

function Base.:(*)(lhs::LadderSumOperator{P, O, S}, rhs::LadderUnitOperator{P, O}) where {P, O, S}
    return LadderSumOperator([t * rhs => a for (t, a) in lhs.terms])
end

function Base.:(*)(lhs::LadderProductOperator{P, O}, rhs::LadderSumOperator{P, O, S}) where {P, O, S}
    return LadderSumOperator([lhs * t => a for (t, a) in rhs.terms])
end

function Base.:(*)(lhs::LadderSumOperator{P, O, S}, rhs::LadderProductOperator{P, O}) where {P, O, S}
    return LadderSumOperator([t * rhs => a for (t, a) in lhs.terms])
end

function Base.:(*)(lhs::LadderSumOperator{P, O, S1}, rhs::LadderSumOperator{P, O, S2}) where {P, O, S1, S2}
    S3 = promote_type(S1, S2)
    terms = Vector{Pair{LadderProductOperator{P, O}, S3}}(undef, length(lhs.terms) * length(rhs.terms))
    i = 0
    for (tl, al) in lhs.terms, (tr, ar) in rhs.terms
        i += 1
        terms[i] = (LadderProductOperator(vcat(tl.factors, tr.factors)) => al * ar)
    end
    return LadderSumOperator(terms)
end


# 3. Addition

function Base.:(+)(lhs::LadderProductOperator{P, O}, rhs::LadderProductOperator{P, O}) where {P, O}
    return LadderSumOperator([lhs=>one(Int), rhs=>one(Int)])
end

function Base.:(-)(lhs::LadderProductOperator{P, O}, rhs::LadderProductOperator{P, O}) where {P, O}
    return LadderSumOperator([lhs=>one(Int), rhs=>-one(Int)])
end

Base.:(+)(lhs::LadderSumOperator{P, O, S1}, rhs::LadderSumOperator{P, O, S2}) where {P, O, S1, S2} = LadderSumOperator(vcat(lhs.terms, rhs.terms))
Base.:(-)(lhs::LadderSumOperator{P, O, S1}, rhs::LadderSumOperator{P, O, S2}) where {P, O, S1, S2} = (lhs) + (-rhs)

function Base.:(+)(lhs::LadderSumOperator{P, O, S}, rhs::LadderProductOperator{P, O}) where {P, O, S}
    return LadderSumOperator([lhs.terms..., rhs=>one(S)])
end

function Base.:(+)(lhs::LadderProductOperator{P, O}, rhs::LadderSumOperator{P, O, S}) where {P, O, S}
    return LadderSumOperator([lhs=>one(S), rhs.terms...])
end

function isequiv(lhs::LadderSumOperator{P, O, S1}, rhs::LadderSumOperator{P, O, S2}) where {P, O, S1, S2}
    return isempty(simplify(lhs - rhs))
end



function Base.adjoint(arg::LadderSumOperator{P, O, S}) where {P, O, S}
    return LadderSumOperator([(adjoint(t), conj(a)) for (t, a) in arg.terms])
end

function LinearAlgebra.ishermitian(arg::LadderSumOperator{P, O, S}) where {P, O, S}
    return isequiv(arg, adjoint(arg))
end
