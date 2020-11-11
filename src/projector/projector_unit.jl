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

    function ParticleProjectorUnitOperator{BR, S}(
        bitmask::Unsigned, bitrow::Unsigned, bitcol::Unsigned,
        parity_bitmask::Unsigned,
        amplitude::Number,
    ) where {BR<:Unsigned, S<:Number}
        if !iszero((~bitmask) & bitrow)
            throw(ArgumentError("bitrow $(string(bitrow, base=2)) has nonzero bit outside bitmask $(string(bitmask, base=2))"))
        elseif !iszero((~bitmask) & bitcol)
            throw(ArgumentError("bitcol $(string(bitcol, base=2)) has nonzero bit outside bitmask $(string(bitmask, base=2))"))
        elseif !iszero(bitmask & parity_bitmask)
            throw(ArgumentError("bitmask and parity_bitmask should have no overlap"))
        end
        return new{BR, S}(bitmask, bitrow, bitcol, parity_bitmask, amplitude)
    end


    function ParticleProjectorUnitOperator{BR, S}(amplitude::Number) where {BR<:Unsigned, S<:Number}
        z = zero(BR)
        return new{BR, S}(z, z, z, z, amplitude)
    end
end

function Base.:(==)(x::ParticleProjectorUnitOperator, y::ParticleProjectorUnitOperator)
    return (
        (x.bitmask == y.bitmask) &&
        (x.bitrow == y.bitrow) &&
        (x.bitcol == y.bitcol) &&
        (x.parity_bitmask == y.parity_bitmask) &&
        (x.amplitude == y.amplitude)
    )
end


function Base.isapprox(
    x::ParticleProjectorUnitOperator, y::ParticleProjectorUnitOperator;
    atol::Real=0, rtol::Real=Base.rtoldefault(x.amplitude,y.amplitude,atol), nans::Bool=false,
)
    return (
        (x.bitmask == y.bitmask) &&
        (x.bitrow == y.bitrow) &&
        (x.bitcol == y.bitcol) &&
        (x.parity_bitmask == y.parity_bitmask) &&
        isapprox(x.amplitude, y.amplitude; atol=atol, rtol=rtol, nans=nans)
    )
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


# Unary

function Base.real(x::ParticleProjectorUnitOperator)
    return ParticleProjectorUnitOperator(x.bitmask, x.bitrow, x.bitcol, x.parity_bitmask, real(x.amplitude))
end

function Base.imag(x::ParticleProjectorUnitOperator)
    return ParticleProjectorUnitOperator(x.bitmask, x.bitrow, x.bitcol, x.parity_bitmask, imag(x.amplitude))
end

function Base.conj(x::ParticleProjectorUnitOperator)
    return ParticleProjectorUnitOperator(x.bitmask, x.bitrow, x.bitcol, x.parity_bitmask, conj(x.amplitude))
end

function Base.adjoint(x::ParticleProjectorUnitOperator)
    return ParticleProjectorUnitOperator(x.bitmask, x.bitcol, x.bitrow, x.parity_bitmask, conj(x.amplitude))
end

function Base.transpose(x::ParticleProjectorUnitOperator)
    return ParticleProjectorUnitOperator(x.bitmask, x.bitcol, x.bitrow, x.parity_bitmask, x.amplitude)
end

function Base.:(-)(x::ParticleProjectorUnitOperator)
    return ParticleProjectorUnitOperator(x.bitmask, x.bitrow, x.bitcol, x.parity_bitmask, -x.amplitude)
end

function Base.:(+)(x::ParticleProjectorUnitOperator)
    return x
end


# Binary (Scale)

function Base.:(*)(lhs::Number, rhs::ParticleProjectorUnitOperator)
    return ParticleProjectorUnitOperator(rhs.bitmask, rhs.bitrow, rhs.bitcol, rhs.parity_bitmask, lhs * rhs.amplitude)
end

function Base.:(*)(lhs::ParticleProjectorUnitOperator, rhs::Number)
    return ParticleProjectorUnitOperator(lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.parity_bitmask, lhs.amplitude * rhs)
end

function Base.:(/)(lhs::ParticleProjectorUnitOperator, rhs::Number)
    return ParticleProjectorUnitOperator(lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.parity_bitmask, lhs.amplitude / rhs)
end

function Base.:(//)(lhs::ParticleProjectorUnitOperator, rhs::Number)
    return ParticleProjectorUnitOperator(lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.parity_bitmask, lhs.amplitude // rhs)
end

function Base.:(\)(lhs::Number, rhs::ParticleProjectorUnitOperator)
    return ParticleProjectorUnitOperator(rhs.bitmask, rhs.bitrow, rhs.bitcol, rhs.parity_bitmask, lhs \ rhs.amplitude)
end


function Base.:(^)(x::ParticleProjectorUnitOperator, y::Integer)
    y < 0 && throw(ArgumentError("power cannot be negative"))
    out = one(x)
    z = x
    while y > 0
        (y & 1 != 0) && (out *= x)
        iszero(out) && return zero(out)
        y >>= 1
        x = x * x
    end
    return out
end
