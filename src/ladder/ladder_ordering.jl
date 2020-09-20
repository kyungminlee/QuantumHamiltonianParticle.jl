struct NormalOrdering{BaseOrdering<:Ordering} <: Ordering end


function Base.lt(
    ::NormalOrdering{BaseOrdering},
    a::LadderUnitOperator{PS, PI, OI},
    b::LadderUnitOperator{PS, PI, OI},
) where {BaseOrdering, PS, PI, OI}
    return isless_normalordering(a,b)
end
