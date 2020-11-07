
Base.:(+)(arg::AbstractParticleLadder) = arg
Base.:(-)(lhs::AbstractParticleLadder{PS, <:Number}, rhs::AbstractParticleLadder{PS, <:Number}) where {PS} = lhs + (-rhs)
Base.:(-)(lhs::AbstractParticleLadder, rhs::Number) = lhs + (-rhs)
Base.:(-)(lhs::Number, rhs::AbstractParticleLadder) = lhs + (-rhs)

# 1. Unary

Base.:(-)(arg::ParticleLadderNull) = arg
Base.:(-)(arg::ParticleLadderUnit) = ParticleLadderSum([ParticleLadderProduct([arg])=>-1])
Base.:(-)(arg::ParticleLadderProduct) = ParticleLadderSum([arg=>-1])
Base.:(-)(arg::ParticleLadderSum) = ParticleLadderSum([t => -a for (t, a) in arg.terms])

# 2. Products

# 2a. Simple
Base.:(*)(lhs::ParticleLadderNull, rhs::Number) = lhs
Base.:(*)(lhs::Number, rhs::ParticleLadderNull) = rhs
Base.:(/)(lhs::ParticleLadderNull, rhs::Number) = lhs
Base.:(\)(lhs::Number, rhs::ParticleLadderNull) = rhs
Base.:(//)(lhs::ParticleLadderNull, rhs::Number) = lhs

Base.:(*)(lhs::ParticleLadderUnit, rhs::Number) = ParticleLadderSum([ParticleLadderProduct([lhs]) => rhs])
Base.:(*)(lhs::Number, rhs::ParticleLadderUnit) = ParticleLadderSum([ParticleLadderProduct([rhs]) => lhs])
Base.:(/)(lhs::ParticleLadderUnit, rhs::Number) = ParticleLadderSum([ParticleLadderProduct([lhs]) => inv(rhs)])
Base.:(\)(lhs::Number, rhs::ParticleLadderUnit) = ParticleLadderSum([ParticleLadderProduct([rhs]) => inv(lhs)])
Base.:(//)(lhs::ParticleLadderUnit, rhs::Number) = ParticleLadderSum([ParticleLadderProduct([lhs]) => 1//rhs])

Base.:(*)(lhs::ParticleLadderProduct, rhs::Number) = ParticleLadderSum([lhs => rhs])
Base.:(*)(lhs::Number, rhs::ParticleLadderProduct) = ParticleLadderSum([rhs => lhs])
Base.:(/)(lhs::ParticleLadderProduct, rhs::Number) = ParticleLadderSum([lhs => inv(rhs)])
Base.:(\)(lhs::Number, rhs::ParticleLadderProduct) = ParticleLadderSum([rhs => inv(lhs)])
Base.:(//)(lhs::ParticleLadderProduct, rhs::Number) = ParticleLadderSum([lhs => 1//rhs])

Base.:(*)(lhs::ParticleLadderSum, rhs::Number) = ParticleLadderSum([(t => a * rhs) for (t, a) in lhs.terms])
Base.:(*)(lhs::Number, rhs::ParticleLadderSum) = ParticleLadderSum([(t => lhs * a) for (t, a) in rhs.terms])
Base.:(/)(lhs::ParticleLadderSum, rhs::Number) = ParticleLadderSum([(t => a / rhs) for (t, a) in lhs.terms])
Base.:(\)(lhs::Number, rhs::ParticleLadderSum) = ParticleLadderSum([(t => lhs \ a) for (t, a) in rhs.terms])
Base.:(//)(lhs::ParticleLadderSum, rhs::Number) = ParticleLadderSum([(t => a // rhs) for (t, a) in lhs.terms])


# 2b. Complex

function Base.:(*)(lhs::ParticleLadderNull{PS}, rhs::AbstractParticleLadder{PS, S}) where {PS, S}
    return lhs
end

function Base.:(*)(lhs::AbstractParticleLadder{PS, S}, rhs::ParticleLadderNull{PS}) where {PS, S}
    return rhs
end


function Base.:(*)(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
    return ParticleLadderProduct([lhs, rhs])
end

function Base.:(*)(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
    return ParticleLadderProduct([lhs, rhs.factors...])
end

function Base.:(*)(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderSum{PS, P, O, S}) where {PS, P, O, S}
    return ParticleLadderSum([lhs * t => a for (t, a) in rhs.terms])
end


function Base.:(*)(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
    return ParticleLadderProduct([lhs.factors..., rhs])
end

function Base.:(*)(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
    return ParticleLadderProduct([lhs.factors..., rhs.factors...])
end

function Base.:(*)(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderSum{PS, P, O, S}) where {PS, P, O, S}
    return ParticleLadderSum([lhs * t => a for (t, a) in rhs.terms])
end


function Base.:(*)(lhs::ParticleLadderSum{PS, P, O, S}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O, S}
    return ParticleLadderSum([t * rhs => a for (t, a) in lhs.terms])
end

function Base.:(*)(lhs::ParticleLadderSum{PS, P, O, S}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O, S}
    return ParticleLadderSum([t * rhs => a for (t, a) in lhs.terms])
end

function Base.:(*)(lhs::ParticleLadderSum{PS, P, O, S1}, rhs::ParticleLadderSum{PS, P, O, S2}) where {PS, P, O, S1, S2}
    S3 = promote_type(S1, S2)
    terms = Vector{Pair{ParticleLadderProduct{PS, P, O}, S3}}(undef, length(lhs.terms) * length(rhs.terms))
    i = 0
    for (tl, al) in lhs.terms, (tr, ar) in rhs.terms
        i += 1
        terms[i] = (ParticleLadderProduct(vcat(tl.factors, tr.factors)) => al * ar)
    end
    return ParticleLadderSum(terms)
end



# 3. Addition

function Base.:(+)(lhs::Number, rhs::ParticleLadderNull)
    return lhs
end

function Base.:(+)(lhs::ParticleLadderNull, rhs::Number)
    return rhs
end


function Base.:(+)(lhs::S, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O, S<:Number}
    return ParticleLadderSum([
        one(ParticleLadderProduct{PS, P, O})=>lhs,
        ParticleLadderProduct([rhs])=>one(S)
    ])
end

function Base.:(+)(lhs::S, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O, S<:Number}
    return ParticleLadderSum([
        one(ParticleLadderProduct{PS, P, O})=>lhs,
        rhs=>one(S),
    ])
end

function Base.:(+)(lhs::Number, rhs::ParticleLadderSum{PS, P, O, S}) where {PS, P, O, S}
    return ParticleLadderSum([
        one(ParticleLadderProduct{PS, P, O})=>lhs,
        rhs.terms...,
    ])
end

function Base.:(+)(lhs::ParticleLadderUnit{PS, P, O}, rhs::S) where {PS, P, O, S<:Number}
    return ParticleLadderSum([
        ParticleLadderProduct([lhs])=>one(S),
        one(ParticleLadderProduct{PS, P, O})=>rhs,
    ])
end

function Base.:(+)(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
    return ParticleLadderSum([
        ParticleLadderProduct([lhs])=>one(Int),
        ParticleLadderProduct([rhs])=>one(Int),
    ])
end

function Base.:(+)(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
    return ParticleLadderSum([
        ParticleLadderProduct([lhs])=>one(Int),
        rhs=>one(Int),
    ])
end

function Base.:(+)(lhs::ParticleLadderUnit{PS, P, O}, rhs::ParticleLadderSum{PS, P, O, S}) where {PS, P, O, S}
    return ParticleLadderSum([
        ParticleLadderProduct([lhs])=>one(Int),
        rhs.terms...
    ])
end

function Base.:(+)(lhs::ParticleLadderProduct{PS, P, O}, rhs::S) where {PS, P, O, S<:Number}
    return ParticleLadderSum([
        lhs=>one(S),
        one(ParticleLadderProduct{PS, P, O})=>rhs,
    ])
end

function Base.:(+)(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
    return ParticleLadderSum([
        lhs=>one(Int),
        ParticleLadderProduct([rhs])=>one(Int),
    ])
end

function Base.:(+)(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
    return ParticleLadderSum([lhs=>one(Int), rhs=>one(Int)])
end

function Base.:(+)(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderSum{PS, P, O, S}) where {PS, P, O, S}
    return ParticleLadderSum([lhs=>one(S), rhs.terms...])
end

function Base.:(+)(lhs::ParticleLadderSum{PS, P, O, S}, rhs::Number) where {PS, P, O, S}
    return ParticleLadderSum([
        lhs.terms...,
        one(ParticleLadderProduct{PS, P, O})=>rhs,
    ])
end

function Base.:(+)(lhs::ParticleLadderSum{PS, P, O, S}, rhs::ParticleLadderUnit{PS, P, O}) where {PS, P, O, S}
    return ParticleLadderSum([
        lhs.terms...,
        ParticleLadderProduct([rhs])=>one(Int),
    ])
end

function Base.:(+)(lhs::ParticleLadderSum{PS, P, O, S}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O, S}
    return ParticleLadderSum([lhs.terms..., rhs=>one(S)])
end

Base.:(+)(lhs::ParticleLadderSum{PS, P, O, S1}, rhs::ParticleLadderSum{PS, P, O, S2}) where {PS, P, O, S1, S2} = ParticleLadderSum(vcat(lhs.terms, rhs.terms))
