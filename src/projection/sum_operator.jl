struct ParticleProjectionSumOperator{BR<:Unsigned, S<:Number}
    terms::Vector{ParticleProjectionUnitOperator{BR, S}}
    function ParticleProjectionSumOperator(terms::AbstractVector{ParticleProjectionUnitOperator{BR, S}}) where {BR, S}
        return new{BR, S}(terms)
    end
end


function Base.promote_rule(
    ::Type{ParticleProjectionSumOperator{BL, SL}},
    ::Type{ParticleProjectionUnitOperator{BR, SR}}
) where {BL, BR, SL, SR}
    B = promote_type(BL, BR)
    S = promote_type(SL, SR)
    return ParticleProjectionSumOperator{B, S}
end

function Base.promote_rule(
    ::Type{ParticleProjectionUnitOperator{BL, SL}},
    ::Type{ParticleProjectionSumOperator{BR, SR}},
) where {BL, BR, SL, SR}
    B = promote_type(BL, BR)
    S = promote_type(SL, SR)
    return ParticleProjectionSumOperator{B, S}
end

function Base.promote_rule(
    ::Type{ParticleProjectionSumOperator{BL, SL}},
    ::Type{ParticleProjectionSumOperator{BR, SR}},
) where {BL, BR, SL, SR}
    B = promote_type(BL, BR)
    S = promote_Type(SL, SR)
    return ParticleProjectionSumOperator{B, S}
end


function Base.convert(
    ::Type{ParticleProjectionSumOperator{B, S}},
    arg::ParticleProjectionUnitOperator) where {B, S}
    projop = ParticleProjectionUnitOperator(
        B(arg.bitmask), B(arg.bitrow), B(arg.bitcol),
        B(arg.parity_bitmask), S(arg.amplitude))
    return ParticleProjectionSumOperator([projop])
end


function Base.zero(::Type{ParticleProjectionSumOperator{BR, S}}) where {BR, S}
    terms = Tuple{ParticleProjectionUnitOperator{BR}, S}[]
    return ParticleProjectionSumOperator(terms)
end


function Base.one(::Type{ParticleProjectionSumOperator{BR, S}}) where {BR, S}
    terms = [one(ParticleProjectionUnitOperator{BR, S})]
    return ParticleProjectionSumOperator(terms)
end

Base.iszero(arg::ParticleProjectionSumOperator) = isempty(arg.terms)


Base.:(*)(x::ParticleProjectionSumOperator, y::Number)  = ParticleProjectionSumOperator([t*y for t in x.terms])
Base.:(/)(x::ParticleProjectionSumOperator, y::Number)  = ParticleProjectionSumOperator([t/y for t in x.terms])
Base.:(//)(x::ParticleProjectionSumOperator, y::Number) = ParticleProjectionSumOperator([t//y for t in x.terms])
Base.:(÷)(x::ParticleProjectionSumOperator, y::Number)  = ParticleProjectionSumOperator([t÷y for t in x.terms])
Base.:(*)(y::Number, x::ParticleProjectionSumOperator)  = ParticleProjectionSumOperator([y*t for t in x.terms])
Base.:(\)(y::Number, x::ParticleProjectionSumOperator)  = ParticleProjectionSumOperator([y\t for t in x.terms])


Base.:(+)(arg::ParticleProjectionSumOperator) = arg
Base.:(-)(arg::ParticleProjectionSumOperator) = ParticleProjectionSumOperator([-t for (t, a) in arg.terms])



Base.:(+)(lhs::ParticleProjectionUnitOperator, rhs::ParticleProjectionUnitOperator) = ParticleProjectionSumOperator([lhs, rhs])
Base.:(+)(lhs::ParticleProjectionUnitOperator, rhs::ParticleProjectionSumOperator) = ParticleProjectionSumOperator([lhs, rhs.terms...])
Base.:(+)(lhs::ParticleProjectionSumOperator, rhs::ParticleProjectionUnitOperator) = ParticleProjectionSumOperator([lhs.terms..., rhs])
Base.:(+)(lhs::ParticleProjectionSumOperator, rhs::ParticleProjectionSumOperator) = ParticleProjectionSumOperator(vcat(lhs.terms, rhs.terms))


function Base.:(*)(x::ParticleProjectionSumOperator{BR, S1}, y::ParticleProjectionSumOperator{BR, S2}) where {BR, S1, S2}
    S3 = promote_type(S1, S2)
    terms = ParticleProjectionUnitOperator{BR, S3}[]
    for t1 in x.terms, t2 in y.terms
        t3 = t1 * t2
        if !iszero(t3)
            push!(terms, t3)
        end
    end
    return ParticleProjectionSumOperator(terms)
end


#=
function Base.:(*)(x::ParticleProjectionSumOperator{BR, S1}, y::ParticleProjectionSumOperator{BR, S2}) where {BR, S1, S2}
    S3 = promote_type(S1, S2)
    terms = Tuple{ParticleProjectionUnitOperator{BR}, S3}[]
    for (t1, a1) in x.terms, (t2, a2) in y.terms
        t3, sgn = projection_product(t1, t2)
        if !iszero(sgn)
        push!(terms, (t3, a1*a2*sgn))
        end
    end
    return ParticleProjectionSumOperator(terms)
end


# TODO: right now it's only one-way
function represent(
    phs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderUnitOperator{ParticleIndex{PS}, <:Integer}
) where {PS, BR, QN}
    iptl = op.particle_index.index
    isite = op.orbital
    ladder = op.ladder

    particle = particle_species(PS, iptl)
    bitmask = get_bitmask(phs, iptl, isite)

    terms = Tuple{ParticleProjectionUnitOperator{BR}, Int}[]
    for n in 1:maxoccupancy(particle)
        lo = BR(n-1) << bitoffset(phs, iptl, isite)
        hi = BR(n) << bitoffset(phs, iptl, isite)
        bitrow, bitcol = (ladder == CREATION) ? (hi, lo) : (lo, hi)

        pmask, prow, pcol, pcheck = zero(BR), zero(BR), zero(BR), zero(BR)
        if isfermion(particle)
            pmask, prow, pcol = bitmask, bitrow, bitcol
            pcheck = get_bitmask(phs, iptl, 1:isite-1)
        end
        push!(terms, (ParticleProjectionUnitOperator(bitmask, bitrow, bitcol, pmask, prow, pcol, pcheck), 1))
    end
    return ParticleProjectionSumOperator(terms)
end


function represent(
    phs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderProductOperator{ParticleIndex{PS}, <:Integer},
) where {PS, BR, QN}
    out = prod(represent(phs, f) for f in op.factors)
    return out
end


function represent(
    phs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderSumOperator{ParticleIndex{PS}, <:Integer, S},
) where {PS, BR, QN, S}
    out = sum(represent(phs, t)*a for (t, a) in op.terms)
    return out
end
=#
