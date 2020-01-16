export normal_order, simplify


function _normal_order(arg::OperatorProduct{<:Tuple{}})
  return (OperatorSum(arg), false)
end


function _normal_order(arg::OperatorProduct{<:Tuple{<:LadderOperator}})
  return (OperatorSum(arg), false)
end


function _normal_order(arg::OperatorProduct{<:Tuple{<:LadderOperator,
                                                    <:LadderOperator,
                                                    Vararg{<:LadderOperator}}})
  factors = arg.factors

  f1, f2 = factors[1], factors[2]
  if maxoccupancy(f1) == 1
    isequal(f1, f2) && return (OperatorSum(), true)
  end

  if isless(f2, f1)
    (tail_op, dirty) = _normal_order(OperatorProduct(f1, factors[3:end]...))
    f2s = OperatorSum(f2) * exchangesign(f1, f2)
    swapped_op = f2s * tail_op
    if isa(f1, AnnihilationOperator) && isa(f2, CreationOperator) && f1.orbital == f2.orbital
      (resid_op, _) = _normal_order(OperatorProduct(factors[3:end]...))
      return (swapped_op + resid_op, true)
    else
      return (swapped_op, true)
    end
  else
    (tail_op, dirty) = _normal_order(OperatorProduct(f2, factors[3:end]...))
    f1s = OperatorSum(f1)
    return (f1s * tail_op, dirty)
  end
end


function _normal_order(arg::OperatorSum{<:Tuple{Vararg{<:Tuple{<:OperatorProduct, <:Number}}}}) ::Tuple{OperatorSum, Bool}
  terms = arg.terms
  out = OperatorSum()
  dirty = false
  for (t, a) in terms
    (t2, d) = _normal_order(t)
    dirty = dirty || d
    out += t2 * a
  end
  return (out, dirty)
end


function normal_order(arg::OperatorSum{<:Tuple{Vararg{<:Tuple{<:OperatorProduct, <:Number}}}}) ::OperatorSum
  isempty(arg.terms) && return arg
  while true
    (arg, dirty) = _normal_order(arg)
    if !dirty
      break
    end
  end
  arg
end


function simplify(arg::OperatorSum{<:Tuple{Vararg{<:Tuple{<:OperatorProduct, <:Number}}}}) ::OperatorSum
  isempty(arg.terms) && return arg
  arg = normal_order(arg)

  terms = sort(Tuple{OperatorProduct, Number}[arg.terms...])
  new_terms = Tuple{OperatorProduct, Number}[]
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
  return OperatorSum(new_terms...)
end


# TODO: rename function
function unify_type(arg::OperatorSum{T}) where T<:Tuple{Vararg{<:Tuple{<:OperatorProduct, <:Number}}}
  isempty(arg.terms) && return arg
  S = promote_type((t.types[2] for t in T.types)...)
  terms = Tuple{OperatorProduct, S}[arg.terms...]
  return OperatorSum(terms...)
end
