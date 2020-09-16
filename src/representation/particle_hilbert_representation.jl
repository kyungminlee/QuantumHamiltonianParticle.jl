struct ParticleHilbertSpaceRepresentation{
    HS<:ParticleHilbertSpace,
}
    function ParticleHilbertSpaceRepresentation(
        hilbert_space::ParticleHilbertSpace{HS, BR, QN}
    )
end
