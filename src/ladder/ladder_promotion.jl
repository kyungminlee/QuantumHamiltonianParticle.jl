# only unit can be converted to product
function Base.promote_rule(::Type{ParticleLadderProduct{PS, P, O}}, ::Type{ParticleLadderUnit{PS, P, O}}) where {PS, P, O}
    return ParticleLadderProduct{PS, P, O}
end

function Base.convert(::Type{ParticleLadderProduct{PS, P, O}}, obj::ParticleLadderUnit{PS, P, O})::ParticleLadderProduct{PS, P, O} where {PS, P, O}
    return ParticleLadderProduct([obj])
end


# any ladder can be converted to sum
function Base.promote_rule(::Type{ParticleLadderSum{PS, P, O, S}}, ::Type{ParticleLadderNull{PS}}) where {PS, P, O, S}
    return ParticleLadderSum{PS, P, O, S}
end

function Base.promote_rule(::Type{ParticleLadderSum{PS, P, O, S}}, ::Type{ParticleLadderUnit{PS, P, O}}) where {PS, P, O, S}
    return ParticleLadderSum{PS, P, O, S}
end

function Base.promote_rule(::Type{ParticleLadderSum{PS, P, O, S}}, ::Type{ParticleLadderProduct{PS, P, O}}) where {PS, P, O, S}
    return ParticleLadderSum{PS, P, O, S}
end

function Base.promote_rule(::Type{ParticleLadderSum{PS, P, O, S1}}, ::Type{ParticleLadderSum{PS, P, O, S2}}) where {PS, P, O, S1, S2}
    S = promote_type(S1, S2)
    return ParticleLadderSum{PS, P, O, S}
end

function Base.convert(::Type{ParticleLadderSum{PS, P, O, S}}, obj::ParticleLadderNull{PS}) where {PS, P, O, S}
    return ParticleLadderSum{PS, P, O, S}()
end

function Base.convert(::Type{ParticleLadderSum{PS, P, O, S}}, obj::ParticleLadderUnit{PS, P, O}) where {PS, P, O, S}
    return ParticleLadderSum([ParticleLadderProduct([obj])=>one(S)])
end

function Base.convert(::Type{ParticleLadderSum{PS, P, O, S}}, obj::ParticleLadderProduct{PS, P, O}) where {PS, P, O, S}
    return ParticleLadderSum([obj=>one(S)])
end

function Base.convert(::Type{ParticleLadderSum{PS, P, O, S}}, obj::ParticleLadderSum{PS, P, O, S2}) where {PS, P, O, S, S2}
    return ParticleLadderSum([t=>convert(S, a) for (t, a) in obj.terms])
end
