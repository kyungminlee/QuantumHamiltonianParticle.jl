import ExactDiagonalization.apply!

function apply!(
    out::SparseState{S1, BR},
    op::AbstractParticleProjectorOperator{BR, S3},
    state::SparseState{S2, BR},
) where {BR, S1, S2, S3}
    for (bvec, ampl) in state.components
        for (newbvec, ampl2) in get_column_iterator(op, bvec)
            a = get(out.components, newbvec, zero(S2))
            out[newbvec] = a + ampl2 * ampl
        end
    end
    return out
end
