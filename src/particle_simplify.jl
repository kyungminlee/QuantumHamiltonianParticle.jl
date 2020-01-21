export normal_order, simplify
export unify_type

import ExactDiagonalization.simplify

function normal_order(arg::OperatorSum{P, O, S}) ::OperatorSum{P, O, S} where {P, O, S}
  isempty(arg.terms) && return arg
  dirty = true
  while dirty
    (arg, dirty) = _normal_order(arg)
  end
  arg
end


function simplify(arg::OperatorSum{P, O, S}) ::OperatorSum{P, O, S} where {P, O, S}
  isempty(arg.terms) && return arg
  arg = normal_order(arg)

  terms = sort(arg.terms)
  new_terms = Tuple{OperatorProduct{P, O}, S}[]
  t, a = Base.first(terms)
  for i in 2:length(terms)
    t2, a2 = terms[i]
    if t2 == t
      a += a2
    else
      if a != 0
        push!(new_terms, (t, a))
      end
      t, a = t2, a2
    end
  end
  if a != 0
    push!(new_terms, (t, a))
  end
  return OperatorSum(new_terms)
end


function _normal_order(arg::OperatorProduct{P, O}) where {P, O}
  if length(arg.factors) <= 1
    return (OperatorSum([(arg, 1)]), false)
  end

  factors = arg.factors
  f1, f2 = factors[1], factors[2]
  if maxoccupancy(f1) == 1 && isequal(f1, f2)
    return (zero(OperatorSum{P,O,Int}), true)
  end

  if isless(f2, f1)
    (tail_op, dirty) = _normal_order(OperatorProduct(vcat(f1, factors[3:end])))
    f2s = OperatorSum([(OperatorProduct([f2]), exchangesign(f1, f2))])
    swapped_op = f2s * tail_op
    if (f1.particle_index == f2.particle_index) && (f1.orbital == f2.orbital) && (f1.ladder == ANNIHILATION) && (f2.ladder == CREATION)
      (resid_op, _) = _normal_order(OperatorProduct(factors[3:end]))
      return (swapped_op + resid_op, true)
    else
      return (swapped_op, true)
    end
  else
    (tail_op, dirty) = _normal_order(OperatorProduct(vcat(f2, factors[3:end])))
    f1s = OperatorSum([(OperatorProduct([f1]), exchangesign(f1, f2))])
    return (f1s * tail_op, dirty)
  end
end


function _normal_order(arg::OperatorSum{P, O, S})::Tuple{OperatorSum{P, O, S}, Bool} where {P, O, S}
  terms = arg.terms
  out = zero(OperatorSum{P, O, S})
  dirty = false
  out_terms = Tuple{OperatorProduct{P, O}, S}[]

  for (t, a) in terms
    (t2, d) = _normal_order(t)
    dirty = dirty || d
    out += t2 * a
  end
  return (out, dirty)
end
