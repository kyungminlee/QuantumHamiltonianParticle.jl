
# 1. Unary

Base.:(+)(arg::LadderUnitOperator{PS, P, O}) where {PS, P, O} = arg
Base.:(-)(arg::LadderUnitOperator{PS, P, O}) where {PS, P, O} = LadderSumOperator([LadderProductOperator([arg])=>-1])

Base.:(+)(arg::LadderProductOperator{PS, P, O}) where {PS, P, O} = arg
Base.:(-)(arg::LadderProductOperator{PS, P, O}) where {PS, P, O} = LadderSumOperator([arg=>-1])

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

function Base.:(*)(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    return LadderProductOperator([lhs, rhs])
end

function Base.:(*)(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return LadderProductOperator([lhs, rhs.factors...])
end

function Base.:(*)(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return LadderSumOperator([lhs * t => a for (t, a) in rhs.terms])
end


function Base.:(*)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    return LadderProductOperator([lhs.factors..., rhs])
end

function Base.:(*)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return LadderProductOperator([lhs.factors..., rhs.factors...])
end

function Base.:(*)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return LadderSumOperator([lhs * t => a for (t, a) in rhs.terms])
end


function Base.:(*)(lhs::LadderSumOperator{PS, P, O, S}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O, S}
    return LadderSumOperator([t * rhs => a for (t, a) in lhs.terms])
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

function Base.:(+)(lhs::S, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O, S<:Number}
    return LadderSumOperator([
        one(LadderProductOperator{PS, P, O})=>lhs,
        LadderProductOperator([rhs])=>one(S)
    ])
end

function Base.:(+)(lhs::S, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O, S<:Number}
    return LadderSumOperator([
        one(LadderProductOperator{PS, P, O})=>lhs,
        rhs=>one(S),
    ])
end

function Base.:(+)(lhs::Number, rhs::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return LadderSumOperator([
        one(LadderProductOperator{PS, P, O})=>lhs,
        rhs.terms...,
    ])
end

function Base.:(+)(lhs::LadderUnitOperator{PS, P, O}, rhs::S) where {PS, P, O, S<:Number}
    return LadderSumOperator([
        LadderProductOperator([lhs])=>one(S),
        one(LadderProductOperator{PS, P, O})=>rhs,
    ])
end

function Base.:(+)(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    return LadderSumOperator([
        LadderProductOperator([lhs])=>one(Int),
        LadderProductOperator([rhs])=>one(Int),
    ])
end

function Base.:(+)(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return LadderSumOperator([
        LadderProductOperator([lhs])=>one(Int),
        rhs=>one(Int),
    ])
end

function Base.:(+)(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return LadderSumOperator([
        LadderProductOperator([lhs])=>one(Int),
        rhs.terms...
    ])
end

function Base.:(+)(lhs::LadderProductOperator{PS, P, O}, rhs::S) where {PS, P, O, S<:Number}
    return LadderSumOperator([
        lhs=>one(S),
        one(LadderProductOperator{PS, P, O})=>rhs,
    ])
end

function Base.:(+)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    return LadderSumOperator([
        lhs=>one(Int),
        LadderProductOperator([rhs])=>one(Int),
    ])
end

function Base.:(+)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return LadderSumOperator([lhs=>one(Int), rhs=>one(Int)])
end

function Base.:(+)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderSumOperator{PS, P, O, S}) where {PS, P, O, S}
    return LadderSumOperator([lhs=>one(S), rhs.terms...])
end

function Base.:(+)(lhs::LadderSumOperator{PS, P, O, S}, rhs::Number) where {PS, P, O, S}
    return LadderSumOperator([
        lhs.terms...,
        one(LadderProductOperator{PS, P, O})=>rhs,
    ])
end

function Base.:(+)(lhs::LadderSumOperator{PS, P, O, S}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O, S}
    return LadderSumOperator([
        lhs.terms...,
        LadderProductOperator([rhs])=>one(Int),
    ])
end

function Base.:(+)(lhs::LadderSumOperator{PS, P, O, S}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O, S}
    return LadderSumOperator([lhs.terms..., rhs=>one(S)])
end

Base.:(+)(lhs::LadderSumOperator{PS, P, O, S1}, rhs::LadderSumOperator{PS, P, O, S2}) where {PS, P, O, S1, S2} = LadderSumOperator(vcat(lhs.terms, rhs.terms))





# function Base.:(-)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
#     return LadderSumOperator([lhs=>one(Int), rhs=>-one(Int)])
# end
