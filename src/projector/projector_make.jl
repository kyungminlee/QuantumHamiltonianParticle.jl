export make_projector_operator

import QuantumHamiltonian.get_space

function make_projector_operator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderUnit{PS, PI, OI},
) where {PS, BR, QN, PI<:Integer, OI<:Integer}
    particle = getspecies(PS, op.particle_index)
    bm  = get_bitmask(hs, op.particle_index, op.orbital)
    pbm = zero(BR)
    if isfermion(particle)
        pbm = get_parity_bitmask(hs, op.particle_index, op.orbital)
        if op.ladder == CREATION
            br = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
            bc = zero(BR)
            return ParticleProjectorUnitOperator(bm, br, bc, pbm, 1)
        else
            br = zero(BR)
            bc = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
            return ParticleProjectorUnitOperator(bm, br, bc, pbm, 1)
        end
    elseif isboson(particle)
        if maxoccupancy(particle) <= 0
            throw(ArgumentError("maximum occupancy cannot be nonpositive")) # COV_EXCL_LINE
        elseif maxoccupancy(particle) == 1
            if op.ladder == CREATION
                br = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
                bc = zero(BR)
                return ParticleProjectorUnitOperator(bm, br, bc, pbm, 1)
            else
                br = zero(BR)
                bc = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
                return ParticleProjectorUnitOperator(bm, br, bc, pbm, 1)
            end
        else
            if op.ladder == CREATION
                let compute_amplitude(cr::Integer, cc::Integer) = sqrt(cr)
                    return ParticleProjectorSumOperator([
                        let cc = cr - 1,
                            br = BR(cr) << bitoffset(hs, op.particle_index, op.orbital),
                            bc = BR(cc) << bitoffset(hs, op.particle_index, op.orbital),
                            ampl = compute_amplitude(cr, cc)
                            ParticleProjectorUnitOperator(bm, br, bc, pbm, ampl)
                        end for cr in 1:maxoccupancy(particle)
                    ])
                end
            else # op.ladder == ANNIHILATION
                let compute_amplitude(cr::Integer, cc::Integer) = sqrt(cc)
                    return ParticleProjectorSumOperator([
                        let cr = cc - 1,
                            br = BR(cr) << bitoffset(hs, op.particle_index, op.orbital),
                            bc = BR(cc) << bitoffset(hs, op.particle_index, op.orbital),
                            ampl = compute_amplitude(cr, cc)
                            ParticleProjectorUnitOperator(bm, br, bc, pbm, ampl)
                        end for cc in 1:maxoccupancy(particle)
                    ])
                end
            end
        end
    elseif isspin(particle)
        M = maxoccupancy(particle)
        if M <= 0
            throw(ArgumentError("maximum occupancy cannot be nonpositive")) # COV_EXCL_LINE
        elseif M == 1
            if op.ladder == CREATION
                br = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
                bc = zero(BR)
                return ParticleProjectorUnitOperator(bm, br, bc, pbm, 1)
            else
                br = zero(BR)
                bc = one(BR) << bitoffset(hs, op.particle_index, op.orbital)
                return ParticleProjectorUnitOperator(bm, br, bc, pbm, 1)
            end
        else
            let compute_amplitude(cr::Integer, cc::Integer) = 0.5*sqrt(2*M*(1+cr+cc) - 4*cr*cc)
                if op.ladder == CREATION
                    return ParticleProjectorSumOperator([
                        let cc = cr - 1,
                            br = BR(cr) << bitoffset(hs, op.particle_index, op.orbital),
                            bc = BR(cc) << bitoffset(hs, op.particle_index, op.orbital),
                            ampl = compute_amplitude(cr, cc)
                            ParticleProjectorUnitOperator(bm, br, bc, pbm, ampl)
                        end for cr in 1:maxoccupancy(particle)
                    ])
                else # op.ladder == ANNIHILATION
                    return ParticleProjectorSumOperator([
                        let cr = cc - 1,
                            br = BR(cr) << bitoffset(hs, op.particle_index, op.orbital),
                            bc = BR(cc) << bitoffset(hs, op.particle_index, op.orbital),
                            ampl = compute_amplitude(cr, cc)
                            ParticleProjectorUnitOperator(bm, br, bc, pbm, ampl)
                        end for cc in 1:maxoccupancy(particle)
                    ])
                end
            end
        end
    else
        throw(ArgumentError("unsupported particle $particle")) # COV_EXCL_LINE
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
    op::ParticleLadderSum{PS, PI, OI, S},
) where {PS, BR, QN, PI<:Integer, OI<:Integer, S}
    return sum(a * make_projector_operator(hs, t) for (t, a) in op.terms)
end

function make_projector_operator(
    hs::HilbertSpaceSector{<:ParticleHilbertSpace, <:Any},
    op...
)
    return make_projector_operator(basespace(hs), op...)
end

# function make_projector_operator(op::ParticleLadderOperatorEmbedding)
#     return make_projector_operator(get_space(op), get_ladder(op))
# end
