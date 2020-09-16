export LadderProductOperator

import LinearAlgebra
import ExactDiagonalization.isequiv

struct LadderProductOperator{PS, P, O}<:AbstractParticleLadderOperator{PS}
    factors::Vector{LadderUnitOperator{PS, P, O}}
    function LadderProductOperator(
        factors::AbstractVector{LadderUnitOperator{PS, P, O}}
    ) where {PS, P, O}
        new{PS, P, O}(factors)
    end
end


function Base.:(==)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return lhs.factors == rhs.factors
end

function Base.:(*)(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    return LadderProductOperator([lhs, rhs])
end

function Base.:(*)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    return LadderProductOperator([lhs.factors..., rhs])
end

function Base.:(*)(lhs::LadderUnitOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return LadderProductOperator([lhs, rhs.factors...])
end

function Base.:(*)(lhs::LadderProductOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return LadderProductOperator([lhs.factors..., rhs.factors...])
end

function Base.convert(::Type{LadderProductOperator{PS, P, O}}, obj::LadderUnitOperator{PS, P, O}) where {PS, P, O}
    return LadderProductOperator([obj])
end

Base.iszero(arg::LadderProductOperator) = false



function Base.isless(lhs::LadderProductOperator{PS, P, O}, rhs::LadderProductOperator{PS, P, O}) where {PS, P, O}
    ll = length(lhs.factors)
    lr = length(rhs.factors)
    return (ll == lr) ? isless(lhs.factors, rhs.factors) : isless(ll, lr)
end

function Base.adjoint(arg::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return LadderProductOperator([adjoint(f) for f in reverse(arg.factors)])
end


function LinearAlgebra.ishermitian(arg::LadderProductOperator{PS, P, O}) where {PS, P, O}
    return isequiv(arg, adjoint(arg))
end
