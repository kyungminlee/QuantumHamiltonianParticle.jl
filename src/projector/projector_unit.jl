export AbstractParticleProjectorOperator
export ParticleProjectorUnitOperator

abstract type AbstractParticleProjectorOperator{BR<:Unsigned, S<:Number}<:AbstractOperator{S} end

struct ParticleProjectorUnitOperator{BR<:Unsigned, S<:Number}<:AbstractParticleProjectorOperator{BR, S}
    bitmask::BR
    bitrow::BR
    bitcol::BR

    parity_bitmask::BR
    amplitude::S

    function ParticleProjectorUnitOperator(bitmask::BR, bitrow::BR, bitcol::BR, parity_bitmask::BR, amplitude::S) where {BR<:Unsigned, S<:Number}
        if !iszero((~bitmask) & bitrow)
            throw(ArgumentError("bitrow $(string(bitrow, base=2)) has nonzero bit outside bitmask $(string(bitmask, base=2))"))
        elseif !iszero((~bitmask) & bitcol)
            throw(ArgumentError("bitcol $(string(bitcol, base=2)) has nonzero bit outside bitmask $(string(bitmask, base=2))"))
        elseif !iszero(bitmask & parity_bitmask)
            throw(ArgumentError("bitmask and parity_bitmask should have no overlap"))
        end
        return new{BR, S}(bitmask, bitrow, bitcol, parity_bitmask, amplitude)
    end
end

function Base.zero(::Type{ParticleProjectorUnitOperator{BR, S}}) where {BR, S}
    z = zero(BR)
    return ParticleProjectorUnitOperator(z, z, z, z, zero(S))
end

function Base.one(::Type{ParticleProjectorUnitOperator{BR, S}}) where {BR, S}
    z = zero(BR)
    return ParticleProjectorUnitOperator(z, z, z, z, one(S))
end

Base.zero(::ParticleProjectorUnitOperator{BR, S}) where {BR, S} = Base.zero(ParticleProjectorUnitOperator{BR, S})
Base.one(::ParticleProjectorUnitOperator{BR, S}) where {BR, S} = Base.one(ParticleProjectorUnitOperator{BR, S})

Base.iszero(arg::ParticleProjectorUnitOperator) = Base.iszero(arg.amplitude)


function Base.:(*)(lhs::ParticleProjectorUnitOperator{BR, S1}, rhs::ParticleProjectorUnitOperator{BR, S2}) where {BR, S1, S2}
    S = promote_type(S1, S2)
    onlylhs_bitmask   =   lhs.bitmask  & (~rhs.bitmask)
    onlyrhs_bitmask   = (~lhs.bitmask) &  rhs.bitmask
    intersect_bitmask =   lhs.bitmask  &  rhs.bitmask
    union_bitmask     =   lhs.bitmask  |  rhs.bitmask

    if (lhs.bitcol & intersect_bitmask) != (rhs.bitrow & intersect_bitmask)
        return zero(ParticleProjectorUnitOperator{BR, S})
    end

    new_bitmask = union_bitmask
    new_bitrow = lhs.bitrow | (rhs.bitrow & onlyrhs_bitmask)
    new_bitcol = (lhs.bitcol & onlylhs_bitmask) | rhs.bitcol

    new_parity_bitmask = (lhs.parity_bitmask ⊻ rhs.parity_bitmask) & (~union_bitmask)

    isparityeven = mod(count_ones((lhs.bitcol & rhs.parity_bitmask) | (rhs.bitrow & lhs.parity_bitmask)), 2) == 0
    new_amplitude = isparityeven ? lhs.amplitude * rhs.amplitude : -lhs.amplitude * rhs.amplitude

    return ParticleProjectorUnitOperator(new_bitmask, new_bitrow, new_bitcol, new_parity_bitmask, new_amplitude)
end

function Base.:(*)(lhs::Number, rhs::ParticleProjectorUnitOperator)
    return ParticleProjectorUnitOperator(rhs.bitmask, rhs.bitrow, rhs.bitcol, rhs.parity_bitmask, lhs * rhs.amplitude)
end

function Base.:(*)(lhs::ParticleProjectorUnitOperator, rhs::Number)
    return ParticleProjectorUnitOperator(lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.parity_bitmask, lhs.amplitude * rhs)
end
