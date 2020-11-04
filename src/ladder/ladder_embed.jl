export ParticleLadderOperatorEmbedding

import ExactDiagonalization.get_row_iterator
import ExactDiagonalization.get_column_iterator
import ExactDiagonalization.get_element

import ExactDiagonalization.get_space
import LatticeTools.embed

struct ParticleLadderOperatorEmbedding{
    HilbertSpaceType<:ParticleHilbertSpace,
    LadderType<:AbstractParticleLadder,
    PS<:ParticleSector,
    S<:Number
} <:AbstractParticleOperator{PS, S}

    hilbert_space::HilbertSpaceType
    operator::LadderType

    function ParticleLadderOperatorEmbedding(
        hilbert_space::ParticleHilbertSpace{PS, BR, QN},
        ladder_operator::AbstractParticleLadder{PS, S},
    ) where {PS, BR, QN, S}
        OP = typeof(ladder_operator)
        return new{ParticleHilbertSpace{PS, BR, QN}, OP, PS, S}(hilbert_space, ladder_operator)
    end
end

get_space(arg::ParticleLadderOperatorEmbedding) = arg.hilbert_space
get_operator(arg::ParticleLadderOperatorEmbedding) = arg.operator

function embed(
    hilbert_space::ParticleHilbertSpace{PS, BR, QN},
    ladder_operator::AbstractParticleLadder{PS, S},
) where {PS, BR, QN, S}
    return ParticleLadderOperatorEmbedding(hilbert_space, ladder_operator)
end

function get_row_iterator(wrap::ParticleLadderOperatorEmbedding, args...; kwargs...)
    return get_row_iterator(wrap.hilbert_space, wrap.operator, args...; kwargs...)
end

function get_column_iterator(wrap::ParticleLadderOperatorEmbedding, args...; kwargs...)
    return get_column_iterator(wrap.hilbert_space, wrap.operator, args...; kwargs...)
end

function get_element(wrap::ParticleLadderOperatorEmbedding, args...; kwargs...)
    return get_element(wrap.hilbert_space, wrap.operator, args...; kwargs...)
end


# Equality
function Base.:(==)(lhs::O, rhs::O) where {O<:ParticleLadderOperatorEmbedding}
    return (lhs.hilbert_space == rhs.hilbert_space) && (lhs.operator == rhs.operator)
end

function ExactDiagonalization.isequiv(
    lhs::ParticleLadderOperatorEmbedding{H, L, P, S1},
    rhs::ParticleLadderOperatorEmbedding{H, L, P, S2},
) where {H, L, P, S1, S2}
    return (get_space(lhs) == get_space(rhs)) &&
        isequiv(get_operator(lhs), get_operator(rhs))
end

function Base.isapprox(
    lhs::ParticleLadderOperatorEmbedding{H, L, P, S1},
    rhs::ParticleLadderOperatorEmbedding{H, L, P, S2};
    atol::Real=0,
    rtol::Real=Base.rtoldefault(S1,S2,atol)
) where {H, L, P, S1, S2}
    return (lhs.hilbert_space == rhs.hilbert_space) &&
        Base.isapprox(lhs.operator, rhs.operator; atol=atol, rtol=rtol)
end


# Unary Operations

Base.iszero(arg::ParticleLadderOperatorEmbedding) = Base.iszero(get_operator(arg))

Base.:(+)(arg::ParticleLadderOperatorEmbedding) = arg
Base.:(-)(arg::ParticleLadderOperatorEmbedding) = ParticleLadderOperatorEmbedding(arg.hilbert_space, -arg.ladder_oerator)
Base.adjoint(arg::ParticleLadderOperatorEmbedding) = ParticleLadderOperatorEmbedding(arg.hilbert_space, Base.adjoint(arg.ladder_operator))
LinearAlgebra.ishermitian(arg::ParticleLadderOperatorEmbedding) = LinearAlgebra.ishermitian(arg.operator)


# Binary Operations

for binop in [:+, :-, :*]
    @eval begin
        function Base.$binop(
            lhs::ParticleLadderOperatorEmbedding{H, L, P, S1},
            rhs::ParticleLadderOperatorEmbedding{H, L, P, S2},
        ) where {H, L, P, S1, S2}
            if lhs.hilbert_space != rhs.hilbert_space
                throw(ArgumentError("lhs and rhs must have the same Hilbert space"))
            end
            return ParticleLadderOperatorEmbedding(lhs.hilbert_space, Base.$binop(lhs.operator, rhs.operator))
        end
    end
end

for binop in [:*, :/, ://, :div]
    @eval begin
        function Base.$binop(lhs::ParticleLadderOperatorEmbedding, rhs::Number)
            return ParticleLadderOperatorEmbedding(lhs.hilbert_space, Base.$binop(lhs.operator, rhs))
        end
    end
end

for binop in [:*, :\]
    @eval begin
        function Base.$binop(lhs::Number, rhs::ParticleLadderOperatorEmbedding)
            return ParticleLadderOperatorEmbedding(rhs.hilbert_space, Base.$binop(lhs, rhs.operator))
        end
    end
end


function apply!(
    out::SparseState{S1, BR},
    op::ParticleLadderOperatorEmbedding{H, L, P, S3},
    state::SparseState{S2, BR},
) where {BR, S1, S2, H, L, P, S3}
    for (bvec, ampl) in state.components
        for (newbvec, ampl_op) in get_column_iterator(get_space(op), get_operator(op), bvec)
            a = get(out.components, newbvec, zero(S2))
            out[newbvec] = a + ampl_op * ampl
        end
    end
    return out
end


function apply!(
    out::SparseState{S1, BR},
    state::SparseState{S2, BR},
    op::ParticleLadderOperatorEmbedding{H, L, P, S3},
) where {BR, S1, S2, H, L, P, S3}

    for (bvec, ampl) in state.components
        for (newbvec, ampl_op) in get_row_iterator(get_space(op), get_operator(op), bvec)
            a = get(out.components, newbvec, zero(S2))
            out[newbvec] = a + ampl * ampl_op
        end
    end

    return out
end


import ExactDiagonalization.simplify

function simplify(arg::ParticleLadderOperatorEmbedding)
    return ParticleLadderOperatorEmbedding(arg.hilbert_space, simplify(arg.operator))
end
