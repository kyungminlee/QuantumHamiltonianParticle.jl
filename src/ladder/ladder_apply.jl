import QuantumHamiltonian.apply!

function apply!(
    out::SparseState{S1, BR},
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::AbstractParticleLadder{PS, S3},
    state::SparseState{S2, BR},
) where {PS, BR, QN, S1, S2, S3}
    for (bvec, ampl) in state.components
        for (newbvec, ampl2) in get_column_iterator(hs, op, bvec)
            a = get(out.components, newbvec, zero(S2))
            out[newbvec] = a + ampl2 * ampl
        end
    end
    return out
end


function apply!(
    out::SparseState{S1, BR},
    op::AbstractParticleLadder{PS, S3},
    state::SparseState{S2, BR},
) where {PS, BR, S1, S2, S3}
    for (bvec, ampl) in state.components
        for (newbvec, ampl2) in get_column_iterator(op, bvec)
            a = get(out.components, newbvec, zero(S2))
            out[newbvec] = a + ampl2 * ampl
        end
    end
    return out
end
