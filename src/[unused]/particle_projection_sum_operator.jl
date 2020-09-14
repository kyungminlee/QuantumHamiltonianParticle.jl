# using ExactDiagonalization

# export AbstractParticleOperator

# export ParticleProjectionSumOperator
# import ExactDiagonalization.represent

# import Base.zero, Base.one
# import Base.*, Base./, Base.\, Base.//, Base.÷
# import Base.+, Base.-

# import Base.promote_rule
# import Base.convert

# struct ParticleProjectionSumOperator{BR<:Unsigned, S<:Number}
#     terms::Vector{Tuple{ParticleProjectionUnitOperator{BR}, S}}
#     function ParticleProjectionSumOperator(terms::AbstractVector{Tuple{ParticleProjectionUnitOperator{BR}, S}}) where {BR, S}
#         return new{BR, S}(terms)
#     end
# end


# function Base.promote_rule(
#     ::Type{ParticleProjectionSumOperator{BL, S}},
#     ::Type{ParticleProjectionUnitOperator{BR}}
# ) where {BL, BR, S}
#     B = promote_type(BL, BR)
#     return ParticleProjectionSumOperator{B, S}
# end

# function Base.promote_rule(
#     ::Type{ParticleProjectionUnitOperator{BL}},
#     ::Type{ParticleProjectionSumOperator{BR, S}},
# ) where {BL, BR, S}
#     B = promote_type(BL, BR)
#     return ParticleProjectionSumOperator{B, S}
# end

# function Base.promote_rule(
#     ::Type{ParticleProjectionSumOperator{BL, SL}},
#     ::Type{ParticleProjectionSumOperator{BR, SR}},
# ) where {BL, BR, SL, SR}
#     B = promote_type(BL, BR)
#     S = promote_Type(SL, SR)
#     return ParticleProjectionSumOperator{B, S}
# end


# function Base.convert(::Type{ParticleProjectionSumOperator{B, S}}, arg::ParticleProjectionUnitOperator) where {B, S}
#     projop = ParticleProjectionUnitOperator(
#         B(arg.bitmask), B(arg.bitrow), B(arg.bitcol),
#         B(arg.pmask), B(arg.prow), B(arg.pcol), B(arg.pcheck))
#     return ParticleProjectionSumOperator([(projop, one(S))])
# end


# function Base.zero(::Type{ParticleProjectionSumOperator{BR, S}}) where {BR, S}
#     terms = Tuple{ParticleProjectionUnitOperator{BR}, S}[]
#     return ParticleProjectionSumOperator(terms)
# end


# function Base.one(::Type{ParticleProjectionSumOperator{BR, S}}) where {BR, S}
#     terms = [(ParticleProjectionUnitOperator{BR}(), one(S))]
#     return ParticleProjectionSumOperator(terms)
# end


# Base.:(*)(x::ParticleProjectionSumOperator, y::Number) = ParticleProjectionSumOperator([(t, a*y) for (t,a) in x.terms])
# Base.:(/)(x::ParticleProjectionSumOperator, y::Number) = ParticleProjectionSumOperator([(t, a/y) for (t,a) in x.terms])
# Base.:(//)(x::ParticleProjectionSumOperator, y::Number) = ParticleProjectionSumOperator([(t, a//y) for (t,a) in x.terms])
# Base.:(÷)(x::ParticleProjectionSumOperator, y::Number) = ParticleProjectionSumOperator([(t, a÷y) for (t,a) in x.terms])
# Base.:(*)(y::Number, x::ParticleProjectionSumOperator) = ParticleProjectionSumOperator([(t, y*a) for (t,a) in x.terms])
# Base.:(\)(y::Number, x::ParticleProjectionSumOperator) = ParticleProjectionSumOperator([(t, y\a) for (t,a) in x.terms])


# Base.:(+)(arg::ParticleProjectionSumOperator) = arg
# Base.:(-)(arg::ParticleProjectionSumOperator) = ParticleProjectionSumOperator([(t, -a) for (t, a) in arg.terms])


# Base.:(+)(lhs::ParticleProjectionSumOperator,  rhs::ParticleProjectionSumOperator) = ParticleProjectionSumOperator(vcat(lhs.terms, rhs.terms))
# Base.:(-)(lhs::ParticleProjectionSumOperator,  rhs::ParticleProjectionSumOperator) = lhs + (-rhs)


# function Base.:(*)(x::ParticleProjectionSumOperator{BR, S1}, y::ParticleProjectionSumOperator{BR, S2}) where {BR, S1, S2}
#     S3 = promote_type(S1, S2)
#     terms = Tuple{ParticleProjectionUnitOperator{BR}, S3}[]
#     for (t1, a1) in x.terms, (t2, a2) in y.terms
#         t3, sgn = projection_product(t1, t2)
#         if !iszero(sgn)
#         push!(terms, (t3, a1*a2*sgn))
#         end
#     end
#     return ParticleProjectionSumOperator(terms)
# end


# # TODO: right now it's only one-way
# function represent(
#     phs::ParticleHilbertSpace{PS, BR, QN},
#     op::LadderUnitOperator{ParticleIndex{PS}, <:Integer}
# ) where {PS, BR, QN}
#     iptl = op.particle_index.index
#     isite = op.orbital
#     ladder = op.ladder

#     particle = particle_species(PS, iptl)
#     bitmask = get_bitmask(phs, iptl, isite)

#     terms = Tuple{ParticleProjectionUnitOperator{BR}, Int}[]
#     for n in 1:maxoccupancy(particle)
#         lo = BR(n-1) << bitoffset(phs, iptl, isite)
#         hi = BR(n) << bitoffset(phs, iptl, isite)
#         bitrow, bitcol = (ladder == CREATION) ? (hi, lo) : (lo, hi)

#         pmask, prow, pcol, pcheck = zero(BR), zero(BR), zero(BR), zero(BR)
#         if isfermion(particle)
#             pmask, prow, pcol = bitmask, bitrow, bitcol
#             pcheck = get_bitmask(phs, iptl, 1:isite-1)
#         end
#         push!(terms, (ParticleProjectionUnitOperator(bitmask, bitrow, bitcol, pmask, prow, pcol, pcheck), 1))
#     end
#     return ParticleProjectionSumOperator(terms)
# end


# function represent(
#     phs::ParticleHilbertSpace{PS, BR, QN},
#     op::LadderProductOperator{ParticleIndex{PS}, <:Integer},
# ) where {PS, BR, QN}
#     out = prod(represent(phs, f) for f in op.factors)
#     return out
# end


# function represent(
#     phs::ParticleHilbertSpace{PS, BR, QN},
#     op::LadderSumOperator{ParticleIndex{PS}, <:Integer, S},
# ) where {PS, BR, QN, S}
#     out = sum(represent(phs, t)*a for (t, a) in op.terms)
#     return out
# end
