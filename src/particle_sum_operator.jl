using ExactDiagonalization

export AbstractParticleOperator

export ParticleSumOperator

import Base.zero, Base.one
import Base.*, Base./, Base.\, Base.//, Base.÷

struct ParticleSumOperator{BR<:Unsigned, S<:Number}
  terms::Vector{Tuple{ParticleProjectionOperator{BR}, S}}
  function ParticleSumOperator(terms::AbstractVector{Tuple{ParticleProjectionOperator{BR}, S}}) where {BR, S}
    return new{BR, S}(terms)
  end
end

function zero(::Type{ParticleSumOperator{BR, S}}) where {BR, S}
  terms = Tuple{ParticleProjectionOperator{BR}, S}[]
  return ParticleSumOperator(terms)
end

function one(::Type{ParticleSumOperator{BR, S}}) where {BR, S}
  terms = [(ParticleProjectionOperator{BR}(), one(S))]
  return ParticleSumOperator(terms)
end

(*)(x::ParticleSumOperator, y::Number) = ParticleSumOperator([(t, a*y) for (t,a) in x.terms])
(/)(x::ParticleSumOperator, y::Number) = ParticleSumOperator([(t, a/y) for (t,a) in x.terms])
(//)(x::ParticleSumOperator, y::Number) = ParticleSumOperator([(t, a//y) for (t,a) in x.terms])
(÷)(x::ParticleSumOperator, y::Number) = ParticleSumOperator([(t, a÷y) for (t,a) in x.terms])

(*)(y::Number, x::ParticleSumOperator) = ParticleSumOperator([(t, y*a) for (t,a) in x.terms])
(\)(y::Number, x::ParticleSumOperator) = ParticleSumOperator([(t, y\a) for (t,a) in x.terms])

import Base.+, Base.-
(+)(arg::ParticleSumOperator) = arg
(-)(arg::ParticleSumOperator) = ParticleSumOperator([(t, -a) for (t, a) in arg.terms])

(+)(lhs::ParticleSumOperator,  rhs::ParticleSumOperator) = ParticleSumOperator(vcat(lhs.terms, rhs.terms))
(-)(lhs::ParticleSumOperator,  rhs::ParticleSumOperator) = lhs + (-rhs)

import Base.promote_rule

function promote_rule(::Type{ParticleSumOperator{BL, S}},
                      ::Type{ParticleProjectionOperator{BR}}) where {BL, BR, S}
  B = promote_type(BL, BR)
  return ParticleSumOperator{B, S}
end

function promote_rule(::Type{ParticleProjectionOperator{BL}},
                      ::Type{ParticleSumOperator{BR, S}}) where {BL, BR, S}
  B = promote_type(BL, BR)
  return ParticleSumOperator{B, S}
end

function promote_rule(::Type{ParticleSumOperator{BL, SL}},
                      ::Type{ParticleSumOperator{BR, SR}}) where {BL, BR, SL, SR}
  B = promote_type(BL, BR)
  S = promote_Type(SL, SR)
  return ParticleSumOperator{B, S}
end


function convert(::Type{ParticleSumOperator{B, S}}, arg::ParticleProjectionOperator) where {B, S}
  projop = ParticleProjectionOperator(
    B(arg.bitmask), B(arg.bitrow), B(arg.bitcol),
    B(arg.pmask), B(arg.prow), B(arg.pcol), B(arg.pcheck))
  return ParticleSumOperator([(projop, one(S))])
end

function (*)(x::ParticleSumOperator{BR, S1}, y::ParticleSumOperator{BR, S2}) where {BR, S1, S2}
  S3 = promote_type(S1, S2)
  terms = Tuple{ParticleProjectionOperator{BR}, S3}[]
  for (t1, a1) in x.terms, (t2, a2) in y.terms
    t3, sgn = projection_product(t1, t2)
    if sgn != 0
      push!(terms, (t3, a1*a2*sgn))
    end
  end
  return ParticleSumOperator(terms)
end


export make_operator
# TODO: right now it's only one-way
function make_operator(phs::ParticleHilbertSpace{PS, BR, QN},
                       op::LadderOperator{ParticleIndex{PS}, Int}) where {PS, BR, QN}
  iptl = op.particle_index.index
  isite = op.orbital
  ladder = op.ladder
  particle = particle_species(PS, iptl)
  bitmask = get_bitmask(phs, iptl, isite)

  terms = Tuple{ParticleProjectionOperator{BR}, Int}[]
  for n in 1:maxoccupancy(particle)
    lo = BR(n-1) << bitoffset(phs, iptl, isite)
    hi = BR(n) << bitoffset(phs, iptl, isite)
    bitrow, bitcol = (if ladder == CREATION
      hi, lo
    else
      lo, hi
    end)
    pmask = zero(BR)
    prow = zero(BR)
    pcol = zero(BR)
    pcheck = zero(BR)

    if isfermion(particle)
      pmask = bitmask
      prow = bitrow
      pcol = bitcol
      pcheck = get_bitmask(phs, iptl, 1:isite-1)
    end
    push!(terms, (ParticleProjectionOperator(bitmask, bitrow, bitcol, pmask, prow, pcol, pcheck), 1))
  end
  return ParticleSumOperator(terms)
end

function make_operator(phs::ParticleHilbertSpace{PS, BR, QN},
                       op::OperatorProduct{ParticleIndex{PS}, Int}) where {PS, BR, QN}
  out = prod(make_operator(phs, f) for f in op.factors)
  return out
end

function make_operator(phs::ParticleHilbertSpace{PS, BR, QN},
                       op::OperatorSum{ParticleIndex{PS}, Int, S}) where {PS, BR, QN, S}
  out = sum(make_operator(phs, t)*a for (t, a) in op.terms)
  return out
end
