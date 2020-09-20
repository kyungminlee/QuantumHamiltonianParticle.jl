# export ParticleProjectorUnitOperator
# struct ParticleProjectorUnitOperator{BR<:Unsigned, S<:Number}
#     bitmask::BR
#     bitrow::BR
#     bitcol::BR

#     # wigner-jordan string
#     pmask::BR
#     prow::BR
#     pcol::BR
#     pcheck::BR

#     amplitude::S

#     function ParticleProjectorUnitOperator(::Type{BR}) where {BR<:Unsigned}
#         z = zero(BR)
#         return new{BR, Int}(z,z,z, z,z,z,z, one(Int))
#     end

#     function ParticleProjectorUnitOperator(
#         bitmask::BR, bitrow::BR, bitcol::BR,
#         pmask::BR, prow::BR, pcol::BR, pcheck::BR, amplitude::S
#     ) where {BR<:Unsigned, S<:Number}
#         return new{BR, S}(bitmask, bitrow, bitcol, pmask, prow, pcol, pcheck, amplitude)
#     end
# end


# function projection_product(
#     x::ParticleProjectorUnitOperator{BR},
#     y::ParticleProjectorUnitOperator{BR},
# )::Tuple{ParticleProjectorUnitOperator{BR}, Int} where {BR}
#     bma = x.bitmask & y.bitmask  # and
#     if (bma & x.bitcol) != (bma & y.bitrow)
#         return (ParticleProjectorUnitOperator{BR}(), 0)
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

#     return (ParticleProjectorUnitOperator(bmask, brow, bcol, pmask, prow, pcol, pcheck), sgn)
# end


# import Base.promote_rule
# function promote_rule(::Type{ParticleProjectorUnitOperator{BL}}, ::Type{ParticleProjectorUnitOperator{BR}}) where {BL, BR}
#     B = promote_type(BL, BR)
#     return ParticleProjectorUnitOperator{B}
# end


# import Base.convert
# function convert(::Type{ParticleProjectorUnitOperator{B}}, arg::ParticleProjectorUnitOperator) where {B}
#     return ParticleProjectorUnitOperator(
#         B(arg.bitmask), B(arg.bitrow), B(arg.bitcol),
#         B(arg.pmask), B(arg.prow), B(arg.pcol), B(arg.pcheck)
#     )
# end
