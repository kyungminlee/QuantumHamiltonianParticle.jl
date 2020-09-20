import ExactDiagonalization.apply!

function apply!(
    out::SparseState{S1, BR},
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::AbstractParticleLadderOperator{PS},
    state::SparseState{S2, BR},
) where {PS, BR, QN, S1, S2}
    for (bvec, ampl) in state.components
        for (newbvec, ampl2) in get_column_iterator(hs, op, bvec)
            a = get(out.components, newbvec, zero(S2))
            out[newbvec] = a + ampl2 * ampl
        end
    end
    return out
end
