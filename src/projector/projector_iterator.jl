import ExactDiagonalization.get_column_iterator
import ExactDiagonalization.get_row_iterator


function get_column_iterator(
    pureop::ParticleProjectorUnitOperator{BR, S},
    bcol::BR2
) where {BR, S, BR2<:Unsigned}
    match::Bool = (bcol & pureop.bitmask) == pureop.bitcol
    element = let
        isparityeven = mod(count_ones(bcol & pureop.parity_bitmask), 2) == 0
        brow = (bcol & ~pureop.bitmask) | pureop.bitrow
        ampl = isparityeven ? pureop.amplitude : -pureop.amplitude
        brow => ampl
    end
    return (element for i in 1:(match ? 1 : 0))
end



function get_column_iterator(sumop::ParticleProjectorSumOperator{BR, S}, bcol::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
    let bcol::BR = BR(bcol),
        terms::Vector{ParticleProjectorUnitOperator{BR, S}} = sumop.terms
        match(pureop::ParticleProjectorUnitOperator{BR, S})::Bool = (bcol & pureop.bitmask) == pureop.bitcol
        function element(pureop::ParticleProjectorUnitOperator{BR, S})::Pair{BR, S}
            isparityeven = mod(count_ones(bcol & pureop.parity_bitmask), 2) == 0
            brow = (bcol & ~pureop.bitmask) | pureop.bitrow
            ampl = isparityeven ? pureop.amplitude : -pureop.amplitude
            return brow => ampl
        end
        return (element(t) for t in terms if match(t))
    end
end
