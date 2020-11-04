export ParticleLadderProduct

import LinearAlgebra
import ExactDiagonalization.isequiv

struct ParticleLadderProduct{PS, P, O}<:AbstractParticleLadder{PS, Int}
    factors::Vector{ParticleLadderUnit{PS, P, O}}
    function ParticleLadderProduct(
        factors::AbstractVector{ParticleLadderUnit{PS, P, O}}
    ) where {PS, P, O}
        new{PS, P, O}(factors)
    end
end


function Base.:(==)(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
    return lhs.factors == rhs.factors
end

function Base.convert(::Type{ParticleLadderProduct{PS, P, O}}, obj::ParticleLadderUnit{PS, P, O}) where {PS, P, O}
    return ParticleLadderProduct([obj])
end

Base.iszero(arg::ParticleLadderProduct) = false
Base.isone(arg::ParticleLadderProduct) = isempty(arg.factors)

function Base.one(::Type{ParticleLadderProduct{PS, PI, OI}}) where {PS, PI, OI}
    return ParticleLadderProduct(ParticleLadderUnit{PS, PI, OI}[])
end
Base.one(::OP) where {OP<:ParticleLadderProduct} = Base.one(OP)

# Orderings


function Base.isless(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
    ll = length(lhs.factors)
    lr = length(rhs.factors)
    return (ll == lr) ? Base.isless(lhs.factors, rhs.factors) : isless(ll, lr)
end

# function isless_localnormalorder(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
#     ll = length(lhs.factors)
#     lr = length(rhs.factors)
#     return (ll == lr) ? isless_localnormalorder(lhs.factors, rhs.factors) : isless(ll, lr)
# end

# function isless_normalorder(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
#     ll = length(lhs.factors)
#     lr = length(rhs.factors)
#     return (ll == lr) ? isless_normalorder(lhs.factors, rhs.factors) : isless(ll, lr)
# end

# function isless_localnormalorder(
#     lhs::Pair{ParticleLadderProduct{PS, P, O}, A},
#     rhs::Pair{ParticleLadderProduct{PS, P, O}, B},
# ) where {PS, P, O, A, B}
#     return isless_localnormalorder(lhs[1], rhs[1])
# end

# function isless_normalorder(lhs::ParticleLadderProduct{PS, P, O}, rhs::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
#     ll = length(lhs.factors)
#     lr = length(rhs.factors)
#     return (ll == lr) ? isless_normalorder(lhs.factors, rhs.factors) : isless(ll, lr)
# end

function Base.adjoint(arg::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
    return ParticleLadderProduct([adjoint(f) for f in reverse(arg.factors)])
end

function LinearAlgebra.ishermitian(arg::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
    return iszero(simplify(arg - adjoint(arg)))
end
