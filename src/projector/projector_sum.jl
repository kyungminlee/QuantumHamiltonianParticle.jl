export ParticleProjectorSumOperator

struct ParticleProjectorSumOperator{BR<:Unsigned, S<:Number}<:AbstractParticleProjectorOperator{BR, S}
    terms::Vector{ParticleProjectorUnitOperator{BR, S}}
    function ParticleProjectorSumOperator(terms::AbstractVector{ParticleProjectorUnitOperator{BR, S}}) where {BR, S}
        return new{BR, S}(terms)
    end
end


function Base.zero(::Type{ParticleProjectorSumOperator{BR, S}}) where {BR, S}
    return ParticleProjectorSumOperator(ParticleProjectorUnitOperator{BR, S}[])
end

function Base.one(::Type{ParticleProjectorSumOperator{BR, S}}) where {BR, S}
    return ParticleProjectorSumOperator([one(ParticleProjectorUnitOperator{BR, S})])
end

Base.iszero(arg::ParticleProjectorSumOperator) = isempty(arg.terms)

Base.:(*)(x::ParticleProjectorSumOperator, y::Number)  = ParticleProjectorSumOperator(x.terms .* y)
Base.:(/)(x::ParticleProjectorSumOperator, y::Number)  = ParticleProjectorSumOperator(x.terms ./ y)
Base.:(//)(x::ParticleProjectorSumOperator, y::Number) = ParticleProjectorSumOperator(x.terms .// y)
Base.:(*)(y::Number, x::ParticleProjectorSumOperator)  = ParticleProjectorSumOperator(y .* x.terms)
Base.:(\)(y::Number, x::ParticleProjectorSumOperator)  = ParticleProjectorSumOperator(y .\ x.terms)

function Base.:(÷)(x::ParticleProjectorSumOperator{BR, S1}, y::S2) where {BR, S1, S2}
    S = promote_type(S1, S2)
    return ParticleProjectorSumOperator([x.terms..., one(ParticleProjectorUnitOperator{BR, S})])
end


Base.:(+)(arg::ParticleProjectorSumOperator) = arg
Base.:(-)(arg::ParticleProjectorSumOperator) = ParticleProjectorSumOperator(-arg.terms)

Base.:(+)(lhs::ParticleProjectorUnitOperator, rhs::ParticleProjectorUnitOperator) = ParticleProjectorSumOperator([lhs, rhs])
Base.:(+)(lhs::ParticleProjectorUnitOperator, rhs::ParticleProjectorSumOperator) = ParticleProjectorSumOperator([lhs, rhs.terms...])
Base.:(+)(lhs::ParticleProjectorSumOperator, rhs::ParticleProjectorUnitOperator) = ParticleProjectorSumOperator([lhs.terms..., rhs])
Base.:(+)(lhs::ParticleProjectorSumOperator, rhs::ParticleProjectorSumOperator) = ParticleProjectorSumOperator(vcat(lhs.terms, rhs.terms))

function Base.:(*)(x::ParticleProjectorSumOperator{BR, S1}, y::ParticleProjectorSumOperator{BR, S2}) where {BR, S1, S2}
    S3 = promote_type(S1, S2)
    terms = ParticleProjectorUnitOperator{BR, S3}[]
    for t1 in x.terms, t2 in y.terms
        t3 = t1 * t2
        !iszero(t3) && push!(terms, t3)
    end
    return ParticleProjectorSumOperator(terms)
end


#=
function Base.:(*)(x::ParticleProjectorSumOperator{BR, S1}, y::ParticleProjectorSumOperator{BR, S2}) where {BR, S1, S2}
    S3 = promote_type(S1, S2)
    terms = Tuple{ParticleProjectorUnitOperator{BR}, S3}[]
    for (t1, a1) in x.terms, (t2, a2) in y.terms
        t3, sgn = projection_product(t1, t2)
        if !iszero(sgn)
        push!(terms, (t3, a1*a2*sgn))
        end
    end
    return ParticleProjectorSumOperator(terms)
end


# TODO: right now it's only one-way
function represent(
    phs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderUnit{ParticleIndex{PS}, <:Integer}
) where {PS, BR, QN}
    iptl = op.particle_index.index
    isite = op.orbital
    ladder = op.ladder

    particle = getspecies(PS, iptl)
    bitmask = get_bitmask(phs, iptl, isite)

    terms = Tuple{ParticleProjectorUnitOperator{BR}, Int}[]
    for n in 1:maxoccupancy(particle)
        lo = BR(n-1) << bitoffset(phs, iptl, isite)
        hi = BR(n) << bitoffset(phs, iptl, isite)
        bitrow, bitcol = (ladder == CREATION) ? (hi, lo) : (lo, hi)

        pmask, prow, pcol, pcheck = zero(BR), zero(BR), zero(BR), zero(BR)
        if isfermion(particle)
            pmask, prow, pcol = bitmask, bitrow, bitcol
            pcheck = get_bitmask(phs, iptl, 1:isite-1)
        end
        push!(terms, (ParticleProjectorUnitOperator(bitmask, bitrow, bitcol, pmask, prow, pcol, pcheck), 1))
    end
    return ParticleProjectorSumOperator(terms)
end


function represent(
    phs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderProduct{ParticleIndex{PS}, <:Integer},
) where {PS, BR, QN}
    out = prod(represent(phs, f) for f in op.factors)
    return out
end


function represent(
    phs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderSum{ParticleIndex{PS}, <:Integer, S},
) where {PS, BR, QN, S}
    out = sum(represent(phs, t)*a for (t, a) in op.terms)
    return out
end
=#
