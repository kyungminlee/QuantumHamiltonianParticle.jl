import Base.==
import Base.+, Base.-, Base.*, Base./, Base.\

export AbstractParticleOperator
export LadderOperator
export CreationOperator, AnnihilationOperator

export particle_species, orbitaltype
export exchangesign
export isfermion
export maxoccupancy

export OperatorSum
export OperatorProduct


abstract type AbstractParticleOperator end

export LadderType, CREATION, ANNIHILATION

@enum LadderType CREATION ANNIHILATION

struct LadderOperator{PI<:ParticleIndex, OI}<:AbstractParticleOperator
  particle_index::PI   # which particle
  orbital::OI          # which orbital
  ladder::LadderType   # creation or annihilation

  function LadderOperator(p::P, o::O, l::LadderType) where {P<:ParticleIndex, O}
    return new{P, O}(p, o, l)
  end
end

function maxoccupancy(arg::LadderOperator{PI, OI}) where {PI, OI}
  return maxoccupancy(particle_species(arg.particle_index))
end

import Base.==
function (==)(lhs::LadderOperator{P, O}, rhs::LadderOperator{P, O}) where {P, O}
  return lhs.particle_index == rhs.particle_index && lhs.orbital == rhs.orbital && lhs.ladder == rhs.ladder
end

import Base.isless
function isless(lhs::LadderOperator{P, O}, rhs::LadderOperator{P, O}) where {P, O}
  lhs.ladder != rhs.ladder && return isless(lhs.ladder, rhs.ladder)

  if lhs.ladder == CREATION # both CREATION
    lhs.particle_index != rhs.particle_index && return isless(lhs.particle_index, rhs.particle_index)
    return isless(lhs.orbital, rhs.orbital)
  else # both ANNIHILATION
    lhs.particle_index != rhs.particle_index && return isless(rhs.particle_index, lhs.particle_index)
    return isless(rhs.orbital, lhs.orbital)
  end
end


function exchangesign(lhs::LadderOperator{P, O}, rhs::LadderOperator{P, O}) where {P, O}
  lhs.particle_index != rhs.particle_index && return 1
  isfermion(particle_species(lhs.particle_index)) && return -1
  return 1
end



struct OperatorProduct{P, O}<:AbstractParticleOperator
  factors::Vector{LadderOperator{P, O}}
  function OperatorProduct(factors::AbstractVector{LadderOperator{P, O}}) where {P, O}
    new{P, O}(factors)
  end
end

function (==)(lhs::OperatorProduct{P, O}, rhs::OperatorProduct{P, O}) where {P, O}
  return lhs.factors == rhs.factors
end

struct OperatorSum{P, O, S<:Number}<:AbstractParticle
  terms::Vector{Tuple{OperatorProduct{P, O}, S}}

  function OperatorSum(op::LadderOperator{P, O}) where {P, O}
    return new{P, O, Int}([(OperatorProduct([op]), 1),])
  end

  function OperatorSum(op::OperatorProduct{P, O}) where {P, O}
    return new{P, O, Int}([(op, 1),])
  end

  function OperatorSum(op::OperatorProduct{P, O}, am::S) where {P, O, S<:Number}
    return new{P, O, S}([(op, am),])
  end

  function OperatorSum(terms::AbstractVector{Tuple{OperatorProduct{P, O}, S}}) where {P, O, S<:Number}
    return new{P, O, S}(terms)
  end
end


function (==)(lhs::OperatorSum{P, O, S}, rhs::OperatorSum{P, O, S}) where {P, O, S}
  return lhs.terms == rhs.terms
end


import Base.one
function one(::Type{OperatorSum{P, O, S}}) where {P, O, S}
  return OperatorSum([(OperatorProduct(LadderOperator{P, O}), one(S))])
end

import Base.zero
function zero(::Type{OperatorSum{P, O, S}}) where {P, O, S}
  return OperatorSum(Tuple{OperatorProduct{P, O}, S}[])
end


import Base.isempty, Base.iszero
isempty(arg::OperatorProduct) = isempty(arg.factors)
isempty(arg::OperatorSum) = isempty(arg.terms)

iszero(arg::OperatorProduct) = false
iszero(arg::OperatorSum) = isempty(arg.terms)


import Base.+, Base.-, Base.*, Base./, Base.\, Base.//, Base.÷

(+)(arg::OperatorSum{P, O, S}) where {P, O, S} = arg
(-)(arg::OperatorSum{P, O, S}) where {P, O, S} = OperatorSum([(t, -a) for (t, a) in arg.terms])

(+)(lhs::OperatorSum{P, O, S1}, rhs::OperatorSum{P, O, S2}) where {P, O, S1, S2} = OperatorSum(vcat(lhs.terms, rhs.terms))
(-)(lhs::OperatorSum{P, O, S1}, rhs::OperatorSum{P, O, S2}) where {P, O, S1, S2} = (lhs) + (-rhs)

(*)(lhs::OperatorSum, rhs::Number) = OperatorSum([(t, a * rhs) for (t, a) in lhs.terms])
(*)(lhs::Number, rhs::OperatorSum) = OperatorSum([(t, lhs * a) for (t, a) in rhs.terms])
(/)(lhs::OperatorSum, rhs::Number) = OperatorSum([(t, a / rhs) for (t, a) in lhs.terms])
(\)(lhs::Number, rhs::OperatorSum) = OperatorSum([(t, lhs \ a) for (t, a) in rhs.terms])
(//)(lhs::OperatorSum, rhs::Number) = OperatorSum([(t, a // rhs) for (t, a) in lhs.terms])
(÷)(lhs::OperatorSum, rhs::Number) = OperatorSum([(t, a ÷ rhs) for (t, a) in lhs.terms])


function isless(lhs::OperatorProduct{P, O}, rhs::OperatorProduct{P, O}) where {P, O}
  ll = length(lhs.factors)
  lr = length(rhs.factors)
  return (ll == lr) ? isless(lhs.factors, rhs.factors) : isless(ll, lr)
end


function isequiv(lhs::OperatorSum{P, O, S1}, rhs::OperatorSum{P, O, S2}) where {P, O, S1, S2}
  return isempty(simplify(lhs - rhs))
end

function (*)(lhs::OperatorSum{P, O, S1}, rhs::OperatorSum{P, O, S2}) where {P, O, S1, S2}
  S3 = promote_type(S1, S2)
  terms = Vector{Tuple{OperatorProduct{P, O}, S3}}(undef, length(lhs.terms) * length(rhs.terms))
  i = 0
  for (tl, al) in lhs.terms, (tr, ar) in rhs.terms
    i += 1
    terms[i] = (OperatorProduct(vcat(tl.factors, tr.factors)), al * ar)
  end
  return OperatorSum(terms)
end



import Base.adjoint
function adjoint(arg::LadderOperator{P, O}) where {P, O}
  new_ladder = arg.ladder == CREATION ? ANNIHILATION : CREATION
  return LadderOperator(arg.particle_index, arg.orbital, new_ladder)
end

function adjoint(arg::OperatorProduct{P, O}) where {P, O}
  return OperatorProduct([adjoint(f) for f in reverse(arg.factors)])
end

function adjoint(arg::OperatorSum{P, O, S}) where {P, O, S}
  return OperatorSum([(adjoint(t), conj(a)) for (t, a) in arg.terms])
end

import LinearAlgebra.ishermitian

ishermitian(arg::LadderOperator) = false
function ishermitian(arg::OperatorProduct{P, O}) where {P, O}
  return isequiv(arg, adjoint(arg))
end

function ishermitian(arg::OperatorSum{P, O, S}) where {P, O, S}
  return isequiv(arg, adjoint(arg))
end

# for fname in [:particle_species, :maxoccupancy, :exchangesign, :isfermion]
#   @eval begin
#     ($fname)(arg::LadderOperator) = ($fname)(arg.particle)
#   end
# end
