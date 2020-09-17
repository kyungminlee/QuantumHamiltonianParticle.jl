export normal_order, simplify
export unify_type

import ExactDiagonalization.simplify


function normal_order(arg::LadderSumOperator{PS, P, O, S})::LadderSumOperator{PS, P, O, S} where {PS, P, O, S}
    isempty(arg.terms) && return arg
    dirty = true
    while dirty
        (arg, dirty) = _normal_order(arg)
    end
    arg
end


function simplify(arg::LadderSumOperator{PS, P, O, S})::LadderSumOperator{PS, P, O, S} where {PS, P, O, S}
    isempty(arg.terms) && return arg
    arg = normal_order(arg)

    terms = sort(arg.terms)
    new_terms = Pair{LadderProductOperator{PS, P, O}, S}[]
    t, a = Base.first(terms)
    for i in 2:length(terms)
        t2, a2 = terms[i]
        if t2 == t
            a += a2
        else
            if a != 0
                push!(new_terms, t => a)
            end
            t, a = t2, a2
        end
    end
    if a != 0
        push!(new_terms, t => a)
    end
    return LadderSumOperator(new_terms)
end


function _normal_order(arg::LadderProductOperator{PS, P, O}) where {PS, P, O}
    if length(arg.factors) <= 1
        return (LadderSumOperator([arg => 1]), false)
    end

    factors = arg.factors
    f1, f2 = factors[1], factors[2]
    if maxoccupancy(f1) == 1 && isequal(f1, f2)
        return (zero(LadderSumOperator{PS, P, O, Int}), true)
    end

    if isless(f2, f1)
        (tail_op, dirty) = _normal_order(LadderProductOperator(vcat(f1, factors[3:end])))
        f2s = LadderSumOperator([LadderProductOperator([f2]) => exchangesign(f1, f2)])
        swapped_op = f2s * tail_op
        if (f1.particle_index == f2.particle_index) && (f1.orbital == f2.orbital) && (f1.ladder == ANNIHILATION) && (f2.ladder == CREATION)
            (resid_op, _) = _normal_order(LadderProductOperator(factors[3:end]))
            return (swapped_op + resid_op, true)
        else
            return (swapped_op, true)
        end
    else
        (tail_op, dirty) = _normal_order(LadderProductOperator(vcat(f2, factors[3:end])))
        f1s = LadderSumOperator([LadderProductOperator([f1]) => 1])
        return (f1s * tail_op, dirty)
    end
end


function _normal_order(arg::LadderSumOperator{PS, P, O, S})::Tuple{LadderSumOperator{PS, P, O, S}, Bool} where {PS, P, O, S}
    out = zero(LadderSumOperator{PS, P, O, S})
    dirty = false
    for (t, a) in arg.terms
        (t2, d) = _normal_order(t)
        dirty = dirty || d
        out += t2 * a
    end
    return (out, dirty)
end
