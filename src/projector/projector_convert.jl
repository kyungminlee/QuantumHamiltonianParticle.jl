function Base.promote_rule(
    ::Type{ParticleProjectorUnitOperator{BL, SL}},
    ::Type{ParticleProjectorUnitOperator{BR, SR}}
) where {BL, BR, SL, SR}
    B = promote_type(BL, BR)
    S = promote_type(SL, SR)
    return ParticleProjectorUnitOperator{B, S}
end


function Base.promote_rule(
    ::Type{ParticleProjectorSumOperator{BL, SL}},
    ::Type{ParticleProjectorUnitOperator{BR, SR}}
) where {BL, BR, SL, SR}
    B = promote_type(BL, BR)
    S = promote_type(SL, SR)
    return ParticleProjectorSumOperator{B, S}
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
    S = promote_Type(SL, SR)
    return ParticleProjectorSumOperator{B, S}
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
    arg::ParticleProjectorUnitOperator) where {B, S}
    projop = ParticleProjectorUnitOperator(
        B(arg.bitmask), B(arg.bitrow), B(arg.bitcol),
        B(arg.parity_bitmask), S(arg.amplitude)
    )
    return ParticleProjectorSumOperator([projop])
end
