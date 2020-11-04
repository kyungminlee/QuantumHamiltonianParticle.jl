struct NormalOrdering{BaseOrdering<:Ordering} <: Ordering end


function Base.lt(
    ::NormalOrdering{BaseOrdering},
    a::ParticleLadderUnit{PS, PI, OI},
    b::ParticleLadderUnit{PS, PI, OI},
) where {BaseOrdering, PS, PI, OI}
    return isless_normalordering(a,b)
end
