# export ParticleProjectionUnitOperator
# struct ParticleProjectionUnitOperator{BR<:Unsigned, S<:Number}
#     bitmask::BR
#     bitrow::BR
#     bitcol::BR

#     # wigner-jordan string
#     pmask::BR
#     prow::BR
#     pcol::BR
#     pcheck::BR

#     amplitude::S

#     function ParticleProjectionUnitOperator(::Type{BR}) where {BR<:Unsigned}
#         z = zero(BR)
#         return new{BR, Int}(z,z,z, z,z,z,z, one(Int))
#     end

#     function ParticleProjectionUnitOperator(
#         bitmask::BR, bitrow::BR, bitcol::BR,
#         pmask::BR, prow::BR, pcol::BR, pcheck::BR, amplitude::S
#     ) where {BR<:Unsigned, S<:Number}
#         return new{BR, S}(bitmask, bitrow, bitcol, pmask, prow, pcol, pcheck, amplitude)
#     end
# end


# function projection_product(
#     x::ParticleProjectionUnitOperator{BR},
#     y::ParticleProjectionUnitOperator{BR},
# )::Tuple{ParticleProjectionUnitOperator{BR}, Int} where {BR}
#     bma = x.bitmask & y.bitmask  # and
#     if (bma & x.bitcol) != (bma & y.bitrow)
#         return (ParticleProjectionUnitOperator{BR}(), 0)
#     end
#     bmo = x.bitmask | y.bitmask  # or
#     bmx = x.bitmask ⊻ bma
#     bmy = y.bitmask ⊻ bma

#     pma = x.pmask & y.pmask
#     pmo = x.pmask | y.pmask
#     pmx = x.pmask ⊻ pma
#     pmy = y.pmask ⊻ pma

#     @assert (pma & x.pcol) == (pma & y.prow)

#     bmask = bmo
#     brow = x.bitrow | (y.bitrow & bmy)
#     bcol = y.bitcol | (x.bitcol & bmx)

#     pmask = pmo
#     prow = x.prow | (y.prow & pmy)
#     pcol = y.pcol | (x.pcol & pmx)

#     parity_sign(bs::BR) = (count_ones(bs) % 2 == 0) ? 1 : -1
#     turn_off(tgt::BR, mask::BR) = tgt ⊻ (tgt & mask)

#     sgn = 1
#     pcheck = x.pcheck
#     sgn *= parity_sign(pcheck & y.prow)
#     pcheck = turn_off(pcheck, y.pmask)
#     pcheck ⊻= y.pcheck
#     sgn *= parity_sign(pcheck & pcol)
#     pcheck = turn_off(pcheck, pmask)

#     return (ParticleProjectionUnitOperator(bmask, brow, bcol, pmask, prow, pcol, pcheck), sgn)
# end


# import Base.promote_rule
# function promote_rule(::Type{ParticleProjectionUnitOperator{BL}}, ::Type{ParticleProjectionUnitOperator{BR}}) where {BL, BR}
#     B = promote_type(BL, BR)
#     return ParticleProjectionUnitOperator{B}
# end


# import Base.convert
# function convert(::Type{ParticleProjectionUnitOperator{B}}, arg::ParticleProjectionUnitOperator) where {B}
#     return ParticleProjectionUnitOperator(
#         B(arg.bitmask), B(arg.bitrow), B(arg.bitcol),
#         B(arg.pmask), B(arg.prow), B(arg.pcol), B(arg.pcheck)
#     )
# end
