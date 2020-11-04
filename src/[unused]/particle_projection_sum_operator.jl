# using ExactDiagonalization

# export AbstractParticleOperator

# export ParticleProjectorSumOperator
# import ExactDiagonalization.represent

# import Base.zero, Base.one
# import Base.*, Base./, Base.\, Base.//, Base.÷
# import Base.+, Base.-

# import Base.promote_rule
# import Base.convert

# struct ParticleProjectorSumOperator{BR<:Unsigned, S<:Number}
#     terms::Vector{Tuple{ParticleProjectorUnitOperator{BR}, S}}
#     function ParticleProjectorSumOperator(terms::AbstractVector{Tuple{ParticleProjectorUnitOperator{BR}, S}}) where {BR, S}
#         return new{BR, S}(terms)
#     end
# end


# function Base.promote_rule(
#     ::Type{ParticleProjectorSumOperator{BL, S}},
#     ::Type{ParticleProjectorUnitOperator{BR}}
# ) where {BL, BR, S}
#     B = promote_type(BL, BR)
#     return ParticleProjectorSumOperator{B, S}
# end

# function Base.promote_rule(
#     ::Type{ParticleProjectorUnitOperator{BL}},
#     ::Type{ParticleProjectorSumOperator{BR, S}},
# ) where {BL, BR, S}
#     B = promote_type(BL, BR)
#     return ParticleProjectorSumOperator{B, S}
# end

# function Base.promote_rule(
#     ::Type{ParticleProjectorSumOperator{BL, SL}},
#     ::Type{ParticleProjectorSumOperator{BR, SR}},
# ) where {BL, BR, SL, SR}
#     B = promote_type(BL, BR)
#     S = promote_Type(SL, SR)
#     return ParticleProjectorSumOperator{B, S}
# end


# function Base.convert(::Type{ParticleProjectorSumOperator{B, S}}, arg::ParticleProjectorUnitOperator) where {B, S}
#     projop = ParticleProjectorUnitOperator(
#         B(arg.bitmask), B(arg.bitrow), B(arg.bitcol),
#         B(arg.pmask), B(arg.prow), B(arg.pcol), B(arg.pcheck))
#     return ParticleProjectorSumOperator([(projop, one(S))])
# end


# function Base.zero(::Type{ParticleProjectorSumOperator{BR, S}}) where {BR, S}
#     terms = Tuple{ParticleProjectorUnitOperator{BR}, S}[]
#     return ParticleProjectorSumOperator(terms)
# end


# function Base.one(::Type{ParticleProjectorSumOperator{BR, S}}) where {BR, S}
#     terms = [(ParticleProjectorUnitOperator{BR}(), one(S))]
#     return ParticleProjectorSumOperator(terms)
# end


# Base.:(*)(x::ParticleProjectorSumOperator, y::Number) = ParticleProjectorSumOperator([(t, a*y) for (t,a) in x.terms])
# Base.:(/)(x::ParticleProjectorSumOperator, y::Number) = ParticleProjectorSumOperator([(t, a/y) for (t,a) in x.terms])
# Base.:(//)(x::ParticleProjectorSumOperator, y::Number) = ParticleProjectorSumOperator([(t, a//y) for (t,a) in x.terms])
# Base.:(÷)(x::ParticleProjectorSumOperator, y::Number) = ParticleProjectorSumOperator([(t, a÷y) for (t,a) in x.terms])
# Base.:(*)(y::Number, x::ParticleProjectorSumOperator) = ParticleProjectorSumOperator([(t, y*a) for (t,a) in x.terms])
# Base.:(\)(y::Number, x::ParticleProjectorSumOperator) = ParticleProjectorSumOperator([(t, y\a) for (t,a) in x.terms])


# Base.:(+)(arg::ParticleProjectorSumOperator) = arg
# Base.:(-)(arg::ParticleProjectorSumOperator) = ParticleProjectorSumOperator([(t, -a) for (t, a) in arg.terms])


# Base.:(+)(lhs::ParticleProjectorSumOperator,  rhs::ParticleProjectorSumOperator) = ParticleProjectorSumOperator(vcat(lhs.terms, rhs.terms))
# Base.:(-)(lhs::ParticleProjectorSumOperator,  rhs::ParticleProjectorSumOperator) = lhs + (-rhs)


# function Base.:(*)(x::ParticleProjectorSumOperator{BR, S1}, y::ParticleProjectorSumOperator{BR, S2}) where {BR, S1, S2}
#     S3 = promote_type(S1, S2)
#     terms = Tuple{ParticleProjectorUnitOperator{BR}, S3}[]
#     for (t1, a1) in x.terms, (t2, a2) in y.terms
#         t3, sgn = projection_product(t1, t2)
#         if !iszero(sgn)
#         push!(terms, (t3, a1*a2*sgn))
#         end
#     end
#     return ParticleProjectorSumOperator(terms)
# end


# # TODO: right now it's only one-way
# function represent(
#     phs::ParticleHilbertSpace{PS, BR, QN},
#     op::ParticleLadderUnit{ParticleIndex{PS}, <:Integer}
# ) where {PS, BR, QN}
#     iptl = op.particle_index.index
#     isite = op.orbital
#     ladder = op.ladder

#     particle = getspecies(PS, iptl)
#     bitmask = get_bitmask(phs, iptl, isite)

#     terms = Tuple{ParticleProjectorUnitOperator{BR}, Int}[]
#     for n in 1:maxoccupancy(particle)
#         lo = BR(n-1) << bitoffset(phs, iptl, isite)
#         hi = BR(n) << bitoffset(phs, iptl, isite)
#         bitrow, bitcol = (ladder == CREATION) ? (hi, lo) : (lo, hi)

#         pmask, prow, pcol, pcheck = zero(BR), zero(BR), zero(BR), zero(BR)
#         if isfermion(particle)
#             pmask, prow, pcol = bitmask, bitrow, bitcol
#             pcheck = get_bitmask(phs, iptl, 1:isite-1)
#         end
#         push!(terms, (ParticleProjectorUnitOperator(bitmask, bitrow, bitcol, pmask, prow, pcol, pcheck), 1))
#     end
#     return ParticleProjectorSumOperator(terms)
# end


# function represent(
#     phs::ParticleHilbertSpace{PS, BR, QN},
#     op::ParticleLadderProduct{ParticleIndex{PS}, <:Integer},
# ) where {PS, BR, QN}
#     out = prod(represent(phs, f) for f in op.factors)
#     return out
# end


# function represent(
#     phs::ParticleHilbertSpace{PS, BR, QN},
#     op::ParticleLadderSum{ParticleIndex{PS}, <:Integer, S},
# ) where {PS, BR, QN, S}
#     out = sum(represent(phs, t)*a for (t, a) in op.terms)
#     return out
# end
