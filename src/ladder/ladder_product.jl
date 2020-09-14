export LadderProductOperator

import LinearAlgebra

struct LadderProductOperator{P, O}<:AbstractParticleOperator
    factors::Vector{LadderUnitOperator{P, O}}
    function LadderProductOperator(factors::AbstractVector{LadderUnitOperator{P, O}}) where {P, O}
        new{P, O}(factors)
    end
end


function Base.:(==)(lhs::LadderProductOperator{P, O}, rhs::LadderProductOperator{P, O}) where {P, O}
    return lhs.factors == rhs.factors
end

function Base.:(*)(lhs::LadderUnitOperator{P, O}, rhs::LadderUnitOperator{P, O}) where {P, O}
    return LadderProductOperator([lhs, rhs])
end

function Base.:(*)(lhs::LadderProductOperator{P, O}, rhs::LadderUnitOperator{P, O}) where {P, O}
    return LadderProductOperator([lhs.factors..., rhs])
end

function Base.:(*)(lhs::LadderUnitOperator{P, O}, rhs::LadderProductOperator{P, O}) where {P, O}
    return LadderProductOperator([lhs, rhs.factors...])
end

function Base.:(*)(lhs::LadderProductOperator{P, O}, rhs::LadderProductOperator{P, O}) where {P, O}
    return LadderProductOperator([lhs.factors..., rhs.factors...])
end

function Base.convert(::Type{LadderProductOperator{P, O}}, obj::LadderUnitOperator{P, O}) where {P, O}
    return LadderProductOperator([obj])
end

Base.iszero(arg::LadderProductOperator) = false



function Base.isless(lhs::LadderProductOperator{P, O}, rhs::LadderProductOperator{P, O}) where {P, O}
    ll = length(lhs.factors)
    lr = length(rhs.factors)
    return (ll == lr) ? isless(lhs.factors, rhs.factors) : isless(ll, lr)
end

function Base.adjoint(arg::LadderProductOperator{P, O}) where {P, O}
    return LadderProductOperator([adjoint(f) for f in reverse(arg.factors)])
end


function LinearAlgebra.ishermitian(arg::LadderProductOperator{P, O}) where {P, O}
    return isequiv(arg, adjoint(arg))
end
