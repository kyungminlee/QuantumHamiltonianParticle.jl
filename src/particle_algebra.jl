import Base.==
import Base.+, Base.-, Base.*, Base./, Base.\

abstract type AbstractOperator end


abstract type LadderOperator{ParticleSpecies<:AbstractParticle, OrbitalType}<:AbstractOperator end


struct CreationOperator{ParticleSpecies<:AbstractParticle, OrbitalType}<:LadderOperator{ParticleSpecies, OrbitalType}
  orbital::OrbitalType
  CreationOperator(::Type{P}, orbital::O) where {P, O} = new{P,O}(orbital)
end


struct AnnihilationOperator{ParticleSpecies<:AbstractParticle, OrbitalType}<:LadderOperator{ParticleSpecies, OrbitalType}
  orbital::OrbitalType
  AnnihilationOperator(::Type{P}, orbital::O) where {P, O} = new{P,O}(orbital)
end


particlespecies(::Type{<:LadderOperator{P, <:Any}}) where P = P
orbitaltype(::Type{<:LadderOperator{<:Any, O}}) where O = O

isfermion(::Type{T}) where T<:LadderOperator = isfermion(particlespecies(T))
isfermion(::T) where T<:LadderOperator = isfermion(particlespecies(T))

maxoccupancy(::Type{T}) where T<:LadderOperator = maxoccupancy(particlespecies(T))
maxoccupancy(::T) where T<:LadderOperator = maxoccupancy(particlespecies(T))

exchangesign(::T1, ::T2) where {T1<:LadderOperator, T2<:LadderOperator} = exchangesign(particlespecies(T1), particlespecies(T2))
exchangesign(::Type{T1}, ::Type{T2}) where {T1<:LadderOperator, T2<:LadderOperator} = exchangesign(particlespecies(T1), particlespecies(T2))



function (==)(lhs::CreationOperator{S, O}, rhs::CreationOperator{S, O}) where {S, O}
  return lhs.orbital == rhs.orbital
end

function (==)(lhs::AnnihilationOperator{S, O}, rhs::AnnihilationOperator{S, O}) where {S, O}
  return lhs.orbital == rhs.orbital
end

struct OperatorProduct{F<:Tuple{Vararg{<:LadderOperator}}}<:AbstractOperator
  factors::F

  # The following constructor is not recommended. Want every operators to satisfy `isbits`:
  # function OperatorProduct{F}(factors) where F
  #   return new{F}(factors)
  # end

  function OperatorProduct(factors::LadderOperator...)
    F = typeof(factors)
    return new{F}(factors)
  end

  function OperatorProduct(factors::Tuple{Vararg{<:LadderOperator}})
    F = typeof(factors)
    return new{F}(factors)
  end
end


struct OperatorSum{T<:Tuple{Vararg{Tuple{<:OperatorProduct, <:Number}}} } <: AbstractOperator
  terms::T

  function OperatorSum()
    return new{Tuple{}}(())
  end

  function OperatorSum(op::LadderOperator)
    terms = ( (OperatorProduct(op), 1), )
    T = typeof(terms)
    return new{T}(terms)
  end

  function OperatorSum(op::OperatorProduct)
    terms = ( (op, 1), )
    T = typeof(terms)
    return new{T}(terms)
  end

  function OperatorSum(op::OperatorProduct, am::Number)
    terms = ( (op, am), )
    T = typeof(terms)
    return new{T}(terms)
  end

  function OperatorSum(op::Vararg{<:Tuple{<:OperatorProduct, <:Number}})
    terms = op
    T = typeof(terms)
    return new{T}(terms)
  end
end


import Base.isempty
isempty(arg::OperatorProduct) = isempty(arg.factors)
isempty(arg::OperatorSum) = isempty(arg.terms)


function prettyprint(arg::CreationOperator)
  print("ψ†")
  print("(")
  print(arg.orbital)
  print(")")
end

function prettyprint(arg::AnnihilationOperator)
  print("ψ")
  print("(")
  print(arg.orbital)
  print(")")
end

function prettyprint(arg::OperatorProduct)
  first = true
  for f in arg.factors
    if !first
      print("⋅")
    end
    prettyprint(f)
    first = false
  end
end

function prettyprint(arg::OperatorSum)
  if isempty(arg)
    print("0")
  else
    t, a = arg.terms[1]
    print("(", a, ")")
    if !isempty(t)
      print("⋅")
      prettyprint(t)
    end
    for (t, a) in arg.terms[2:end]
      print(" + (", a, ")")
      if !isempty(t)
        print("⋅")
        prettyprint(t)
      end
    end
  end
end


(+)(arg::OperatorSum) = arg
(-)(arg::OperatorSum) = OperatorSum(((t, -a) for (t, a) in arg.terms)...)

(+)(lhs::OperatorSum, rhs::OperatorSum) = OperatorSum(lhs.terms..., rhs.terms...)

(-)(lhs::OperatorSum, rhs::OperatorSum) = (lhs) + (-rhs)

(*)(lhs::OperatorSum, rhs::Number) = OperatorSum(((t, a * rhs) for (t, a) in lhs.terms)...)
(*)(lhs::Number, rhs::OperatorSum) = OperatorSum(((t, lhs * a) for (t, a) in rhs.terms)...)
(/)(lhs::OperatorSum, rhs::Number) = OperatorSum(((t, a / rhs) for (t, a) in lhs.terms)...)
(\)(lhs::Number, rhs::OperatorSum) = OperatorSum(((t, lhs \ a) for (t, a) in rhs.terms)...)



function (*)(lhs::OperatorSum, rhs::OperatorSum)
  terms = Vector{Any}(undef, length(lhs.terms) * length(rhs.terms))
  i = 0
  for (tl, al) in lhs.terms, (tr, ar) in rhs.terms
    i += 1
    terms[i] = (OperatorProduct(tl.factors..., tr.factors...), al * ar)
  end
  return OperatorSum(terms...)
end


import Base.isless

isless(lhs::CreationOperator, rhs::AnnihilationOperator) = true
isless(lhs::AnnihilationOperator, rhs::CreationOperator) = false

function isless(lhs::CreationOperator{S1, O}, rhs::CreationOperator{S2, O}) where {S1, S2, O}
  if S1 === S2
    return isless(lhs.orbital, rhs.orbital)
  else
    return isless(speciessymb(S1), speciessymb(S2))
  end
end

function isless(lhs::AnnihilationOperator{S1, O}, rhs::AnnihilationOperator{S2, O}) where {S1, S2, O}
  if S1 === S2
    return isless(rhs.orbital, lhs.orbital)
  else
    return isless(speciessymb(S2), speciessymb(S1))
  end
end

function isless(lhs::OperatorProduct, rhs::OperatorProduct) where Index
  ll = length(lhs.factors)
  lr = length(rhs.factors)
  return (ll == lr) ? isless(lhs.factors, rhs.factors) : isless(ll, lr)
end


function isequiv(lhs::OperatorSum, rhs::OperatorSum) where {Index, S1, S2}
  return isempty(simplify(lhs - rhs))
end


