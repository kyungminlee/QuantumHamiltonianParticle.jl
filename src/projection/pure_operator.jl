export AbstractParticleProjectionOperator
export ParticleProjectionUnitOperator

abstract type AbstractParticleProjectionOperator end

Base.:(-)(lhs::AbstractParticleProjectionOperator, rhs::AbstractParticleProjectionOperator) = lhs + (-rhs)

struct ParticleProjectionUnitOperator{BR<:Unsigned, S<:Number} <: AbstractParticleProjectionOperator
    bitmask::BR
    bitrow::BR
    bitcol::BR

    parity_bitmask::BR
    amplitude::S

    function ParticleProjectionUnitOperator(bitmask::BR, bitrow::BR, bitcol::BR, parity_bitmask::BR, amplitude::S) where {BR<:Unsigned, S<:Number}
        if bitmask & parity_bitmask != 0
            throw(ArgumentError("bitmask and parity_bitmask should have no overlap"))
        end
        return new{BR, S}(bitmask, bitrow, bitcol, parity_bitmask, amplitude)
    end
end

function Base.zero(::Type{ParticleProjectionUnitOperator{BR, S}}) where {BR, S}
    z = zero(BR)
    return ParticleProjectionUnitOperator(z, z, z, z, zero(S))
end

function Base.one(::Type{ParticleProjectionUnitOperator{BR, S}}) where {BR, S}
    z = zero(BR)
    return ParticleProjectionUnitOperator(z, z, z, z, one(S))
end

Base.iszero(arg::ParticleProjectionUnitOperator) = Base.iszero(arg.amplitude)


function Base.:(*)(lhs::ParticleProjectionUnitOperator{BR, S1}, rhs::ParticleProjectionUnitOperator{BR, S2}) where {BR, S1, S2}
    S = promote_type(S1, S2)
    onlylhs_bitmask   =   lhs.bitmask  & (~rhs.bitmask)
    onlyrhs_bitmask   = (~lhs.bitmask) &  rhs.bitmask
    intersect_bitmask =   lhs.bitmask  &  rhs.bitmask
    union_bitmask     =   lhs.bitmask  |  rhs.bitmask

    if (lhs.bitcol & intersect_bitmask) != (rhs.bitrow & intersect_bitmask)
        return zero(ParticleProjectionUnitOperator{BR, S})
    end

    new_bitmask = union_bitmask
    new_bitrow = lhs.bitrow | (rhs.bitrow & onlyrhs_bitmask)
    new_bitcol = (lhs.bitcol & onlylhs_bitmask) | rhs.bitcol

    new_parity_bitmask = (lhs.parity_bitmask ⊻ rhs.parity_bitmask) & (~union_bitmask)

    isparityeven = mod(count_ones((lhs.bitcol & rhs.parity_bitmask) | (rhs.bitrow & lhs.parity_bitmask)), 2) == 0
    new_amplitude = isparityeven ? lhs.amplitude * rhs.amplitude : -lhs.amplitude * rhs.amplitude

    return ParticleProjectionUnitOperator(new_bitmask, new_bitrow, new_bitcol, new_parity_bitmask, new_amplitude)
end

function Base.:(*)(lhs::Number, rhs::ParticleProjectionUnitOperator)
    return ParticleProjectionUnitOperator(rhs.bitmask, rhs.bitrow, rhs.bitcol, rhs.parity_bitmask, lhs * rhs.amplitude)
end

function Base.:(*)(lhs::ParticleProjectionUnitOperator, rhs::Number)
    return ParticleProjectionUnitOperator(lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.parity_bitmask, lhs.amplitude * rhs)
end


export make_projection_operator
function make_projection_operator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderUnitOperator{<:Integer, <:Integer},
) where {PS, BR, QN}
    particle = particle_species(PS, op.particle_index)
    bm  = get_bitmask(hs, op.particle_index, op.orbital)
    if isfermion(particle)
        pbm = get_parity_mask(hs, op.particle_index, op.orbital)
        if op.ladder == CREATION
            br = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
            bc = zero(BR)
            return ParticleProjectionUnitOperator(bm, br, bc, pbm, 1.0)
        else
            br = zero(BR)
            bc = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
            return ParticleProjectionUnitOperator(bm, br, bc, pbm, 1.0)
        end
    elseif isboson(particle)
        if maxoccupancy(particle) <= 0
            return NullOperator()
        elseif maxoccupancy(particle) == 1
            pbm = zero(BR)
            if op.ladder == CREATION
                br = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
                bc = zero(BR)
                return ParticleProjectionUnitOperator(bm, br, bc, pbm, 1.0)
            else
                br = zero(BR)
                bc = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
                return ParticleProjectionUnitOperator(bm, br, bc, pbm, 1.0)
            end
        else
            @error "make_projection_operator for bosons and other particles not implemented yet"
        end
    else
        @error "unsupported particle $particle"
    end
end

function make_projection_operator(
    hs::ParticleHilbertSpace,
    op::LadderProductOperator{<:Integer, <:Integer},
)
    return prod(make_projection_operator(hs, f) for f in op.factors)
end

function make_projection_operator(
    hs::ParticleHilbertSpace,
    op::LadderSumOperator{<:Integer, <:Integer, S},
) where {S}
    return sum(a * make_projection_operator(hs, t) for (t, a) in op.terms)
end
