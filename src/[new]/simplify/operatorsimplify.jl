

# Convention:
# wave function
#   c†(1) c†(2) c†(3) |0>
#   Note: Wigner-Jordan:
#
#     111111-zz            |0 0| |1  0| = |0 0|
#     11-zzzzzz            |1 0| |0 -1|   |1 0|
#     ---------
#     11-zzz-11
#
# operator
#  c†(1) c†(2) c†(3)  c(3) c(2) c(1)
#
function isless(lhs::ParticleCreationOperator{P}, rhs::ParticleCreationOperator{P}) where P
    if lhs.particleindex < rhs.particleindex
        return true
    elseif lhs.particleindex > rhs.particleindex
        return false
    else
        return lhs.siteindex < rhs.siteindex
    end
end

function isless(lhs::ParticleAnnihilationOperator{P}, rhs::ParticleAnnihilationOperator{P}) where P
    if lhs.particleindex < rhs.particleindex
        return false
    elseif lhs.particleindex > rhs.particleindex
        return true
    else
        return lhs.siteindex > rhs.siteindex
    end
end


simplify(op::ParticleCreationOperator) = op
simplify(op::ParticleAnnihilationOperator) = op

function simplify(op::AbstractParticleOperator)
    @warn "simplify not yet implemented for $(typeof(op))"
    return op
end
