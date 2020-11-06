export normal_order, simplify
export unify_type

import ExactDiagonalization.simplify

# struct LocalNormalOrdering <: Ordering end
# struct NormalOrdering <: Ordering end

# function lt(
#     ::LocalNormalOrdering,
#     a::ParticleLadderUnit{PS, PI, OI},
#     b::ParticleLadderUnit{PS, PI, OI})
#     isless_localnormalorder(a, b)
# end


function normal_order(arg::ParticleLadderSum{PS, P, O, S})::ParticleLadderSum{PS, P, O, S} where {PS, P, O, S}
    isempty(arg.terms) && return arg
    dirty = true
    while dirty
        (arg, dirty) = _normal_order(arg)
    end
    arg
end

function simplify(arg::ParticleLadderProduct{PS, P, O})::ParticleLadderSum{PS, P, O, Int} where {PS, P, O}
    return simplify(ParticleLadderSum([arg=>1]))
end

simplify(arg::ParticleLadderUnit) = arg

function simplify(arg::ParticleLadderSum{PS, P, O, S})::ParticleLadderSum{PS, P, O, S} where {PS, P, O, S}
    isempty(arg.terms) && return arg
    arg = normal_order(arg)
    isempty(arg.terms) && return arg

    terms = sort(arg.terms)
    new_terms = Pair{ParticleLadderProduct{PS, P, O}, S}[]
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
    return ParticleLadderSum(new_terms)
end


function _normal_order(arg::ParticleLadderProduct{PS, P, O}) where {PS, P, O}
    if length(arg.factors) <= 1
        return (ParticleLadderSum([arg => 1]), false)
    end

    factors = arg.factors
    f1, f2 = factors[1], factors[2]
    if maxoccupancy(f1) == 1 && isequal(f1, f2)
        return (zero(ParticleLadderSum{PS, P, O, Int}), true)
    end

    if isless(f2, f1)
        (tail_op, dirty) = _normal_order(ParticleLadderProduct(vcat(f1, factors[3:end])))
        f2s = ParticleLadderSum([ParticleLadderProduct([f2]) => exchangesign(f1, f2)])
        swapped_op = f2s * tail_op
        if (f1.particle_index == f2.particle_index) && (f1.orbital == f2.orbital) && (f1.ladder == ANNIHILATION) && (f2.ladder == CREATION)
            (resid_op, _) = _normal_order(ParticleLadderProduct(factors[3:end]))
            return (swapped_op + resid_op, true)
        else
            return (swapped_op, true)
        end
    else
        (tail_op, dirty) = _normal_order(ParticleLadderProduct(vcat(f2, factors[3:end])))
        f1s = ParticleLadderSum([ParticleLadderProduct([f1]) => 1])
        return (f1s * tail_op, dirty)
    end
end


function _normal_order(arg::ParticleLadderSum{PS, P, O, S})::Tuple{ParticleLadderSum{PS, P, O, S}, Bool} where {PS, P, O, S}
    out = zero(ParticleLadderSum{PS, P, O, S})
    dirty = false
    for (t, a) in arg.terms
        (t2, d) = _normal_order(t)
        dirty = dirty || d
        out += t2 * a
    end
    return (out, dirty)
end
