export ParticleProjectionOperator
struct ParticleProjectionOperator{BR<:Unsigned}
  bitmask::BR
  bitrow::BR
  bitcol::BR

  # wigner-jordan string
  pmask::BR
  prow::BR
  pcol::BR
  pcheck::BR
  function ParticleProjectionOperator(binary_type::Type{BR}) where {BR<:Unsigned}
    z = zero(BR)
    return new{BR}(z,z,z, z,z,z,z)
  end

  function ParticleProjectionOperator(bitmask::BR, bitrow::BR, bitcol::BR,
                                      pmask::BR, prow::BR, pcol::BR, pcheck::BR) where {BR<:Unsigned}
    return new{BR}(bitmask, bitrow, bitcol, pmask, prow, pcol, pcheck)
  end
end

function projection_product(x::ParticleProjectionOperator{BR},
                            y::ParticleProjectionOperator{BR})::Tuple{ParticleProjectionOperator{BR}, Int} where {BR}
  bma = x.bitmask & y.bitmask  # and
  if (bma & x.bitcol) != (bma & y.bitrow)
    return (ParticleProjectionOperator{BR}(), 0)
  end
  bmo = x.bitmask | y.bitmask  # or
  bmx = x.bitmask ⊻ bma
  bmy = y.bitmask ⊻ bma

  pma = x.pmask & y.pmask
  pmo = x.pmask | y.pmask
  pmx = x.pmask ⊻ pma
  pmy = y.pmask ⊻ pma

  @assert (pma & x.pcol) == (pma & y.prow)

  bmask = bmo;
  brow = x.bitrow | (y.bitrow & bmy);
  bcol = y.bitcol | (x.bitcol & bmx);

  pmask = pmo;
  prow = x.prow | (y.prow & pmy);
  pcol = y.pcol | (x.pcol & pmx);

  parity_sign(bs::BR) = (count_ones(bs) % 2 == 0) ? 1 : -1
  turn_off(tgt::BR, mask::BR) = tgt ⊻ (tgt & mask)

  sgn = 1
  pcheck = x.pcheck
  sgn *= parity_sign(pcheck & y.prow)
  pcheck = turn_off(pcheck, y.pmask)
  pcheck ⊻= y.pcheck
  sgn *= parity_sign(pcheck & pcol)
  pcheck = turn_off(pcheck, pmask)

  return (ParticleProjectionOperator(bmask, brow, bcol, pmask, prow, pcol, pcheck), sgn)
end


import Base.promote_rule
function promote_rule(::Type{ParticleProjectionOperator{BL}},
                      ::Type{ParticleProjectionOperator{BR}}) where {BL, BR}
  B = promote_type(BL, BR)
  return ParticleProjectionOperator{B}
end

import Base.convert
function convert(::Type{ParticleProjectionOperator{B}}, arg::ParticleProjectionOperator) where {B}
  return ParticleProjectionOperator(
    B(arg.bitmask), B(arg.bitrow), B(arg.bitcol),
    B(arg.pmask), B(arg.prow), B(arg.pcol), B(arg.pcheck)
  )
end
