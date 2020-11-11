function Base.promote_rule(
    ::Type{ParticleProjectorUnitOperator{B, SL}},
    ::Type{SR}
) where {B, SL, SR<:Number}
    S = promote_type(SL, SR)
    return ParticleProjectorUnitOperator{B, S}
end

function Base.promote_rule(
    ::Type{ParticleProjectorUnitOperator{B, S}},
    ::Type{<:NullOperator}
) where {B, S}
    return ParticleProjectorUnitOperator{B, S}
end



function Base.promote_rule(
    ::Type{ParticleProjectorUnitOperator{BL, SL}},
    ::Type{ParticleProjectorUnitOperator{BR, SR}}
) where {BL, BR, SL, SR}
    B = promote_type(BL, BR)
    S = promote_type(SL, SR)
    return ParticleProjectorUnitOperator{B, S}
end


function Base.promote_rule(
    ::Type{ParticleProjectorUnitOperator{BL, SL}},
    ::Type{ParticleProjectorSumOperator{BR, SR}},
) where {BL, BR, SL, SR}
    B = promote_type(BL, BR)
    S = promote_type(SL, SR)
    return ParticleProjectorSumOperator{B, S}
end


function Base.promote_rule(
    ::Type{ParticleProjectorSumOperator{BL, SL}},
    ::Type{ParticleProjectorSumOperator{BR, SR}},
) where {BL, BR, SL, SR}
    B = promote_type(BL, BR)
    S = promote_type(SL, SR)
    return ParticleProjectorSumOperator{B, S}
end


function Base.convert(
    ::Type{ParticleProjectorUnitOperator{B, SL}},
    arg::SR,
) where {B, SL, SR}
    S = promote_type(SL, SR)
    z = zero(B)
    return ParticleProjectorUnitOperator(z, z, z, z, S(arg))
end


function Base.convert(
    ::Type{ParticleProjectorUnitOperator{BR, S}},
    ::NullOperator
) where {BR, S}
    z = zero(BR)
    return ParticleProjectorUnitOperator(z, z, z, z, zero(S))
end


function Base.convert(
    ::Type{ParticleProjectorUnitOperator{B, S}},
    arg::ParticleProjectorUnitOperator,
) where {B, S}
    return ParticleProjectorUnitOperator(
        B(arg.bitmask), B(arg.bitrow), B(arg.bitcol),
        B(arg.parity_bitmask), S(arg.amplitude)
    )
end


function Base.convert(
    ::Type{ParticleProjectorSumOperator{B, S}},
    arg::ParticleProjectorUnitOperator,
) where {B, S}
    projop = ParticleProjectorUnitOperator(
        B(arg.bitmask), B(arg.bitrow), B(arg.bitcol),
        B(arg.parity_bitmask), S(arg.amplitude)
    )
    return ParticleProjectorSumOperator([projop])
end



function Base.convert(
    ::Type{ParticleProjectorSumOperator{B, S}},
    arg::ParticleProjectorSumOperator,
) where {B, S}
    terms = ParticleProjectorUnitOperator{B, S}[]
    append!(terms, arg.terms)
    return ParticleProjectorSumOperator(terms)
end
