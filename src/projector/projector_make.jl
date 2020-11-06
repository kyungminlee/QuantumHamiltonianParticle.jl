export make_projector_operator

import ExactDiagonalization.get_space

function make_projector_operator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderUnit{PS, PI, OI},
) where {PS, BR, QN, PI<:Integer, OI<:Integer}
    particle = getspecies(PS, op.particle_index)
    bm  = get_bitmask(hs, op.particle_index, op.orbital)
    if isfermion(particle)
        pbm = get_parity_bitmask(hs, op.particle_index, op.orbital)
        if op.ladder == CREATION
            br = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
            bc = zero(BR)
            return ParticleProjectorUnitOperator(bm, br, bc, pbm, 1.0)
        else
            br = zero(BR)
            bc = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
            return ParticleProjectorUnitOperator(bm, br, bc, pbm, 1.0)
        end
    elseif isboson(particle)
        if maxoccupancy(particle) <= 0
            return NullOperator()
        elseif maxoccupancy(particle) == 1
            pbm = zero(BR)
            if op.ladder == CREATION
                br = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
                bc = zero(BR)
                return ParticleProjectorUnitOperator(bm, br, bc, pbm, 1.0)
            else
                br = zero(BR)
                bc = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
                return ParticleProjectorUnitOperator(bm, br, bc, pbm, 1.0)
            end
        else
            @error "make_projector_operator for bosons and other particles not implemented yet"
            # TODO: Implement ParticleProjectorSumOperator
        end
    else
        @error "unsupported particle $particle"
    end
end

function make_projector_operator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderProduct{PS, PI, OI},
) where {PS, BR, QN, PI<:Integer, OI<:Integer}
    return prod(make_projector_operator(hs, f) for f in op.factors)
end

function make_projector_operator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderSum{PS, PI, OI, S}
) where {PS, BR, QN, PI<:Integer, OI<:Integer, S}
    return sum(a * make_projector_operator(hs, t) for (t, a) in op.terms)
end

function make_projector_operator(
    hs::HilbertSpaceSector{<:ParticleHilbertSpace, <:Any},
    op...
)
    return make_projector_operator(basespace(hs), op...)
end

function make_projector_operator(op::ParticleLadderOperatorEmbedding)
    return make_projector_operator(get_space(op), get_ladder(op))
end
