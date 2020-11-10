import ExactDiagonalization.simplify


function Base.isless(
    lhs::ParticleProjectorUnitOperator{BR, S1},
    rhs::ParticleProjectorUnitOperator{BR, S2},
) where {BR, S1, S2}
    lhs.bitmask < rhs.bitmask && return true
    lhs.bitmask > rhs.bitmask && return false
    lhs.bitrow < rhs.bitrow && return true
    lhs.bitrow > rhs.bitrow && return false
    lhs.bitcol < rhs.bitcol && return true
    lhs.bitcol > rhs.bitcol && return false
    lhs.parity_bitmask < rhs.parity_bitmask && return true
    lhs.parity_bitmask > rhs.parity_bitmask && return false
    real(lhs.amplitude) < real(rhs.amplitude) && return true
    real(lhs.amplitude) > real(rhs.amplitude) && return false
    imag(lhs.amplitude) < imag(rhs.amplitude) && return true
    imag(lhs.amplitude) > imag(rhs.amplitude) && return false
    return false
end


function ExactDiagonalization.simplify(
    x::ParticleProjectorUnitOperator{BR, S},
    tol::Real=Base.rtoldefault(real(S))
) where {BR, S}
    if isapprox(x.amplitude, zero(S); atol=tol)
        return NullOperator()
    else
        if S <: Complex && isapprox(imag(x.amplitude), zero(real(S)); atol=tol)
            return real(x)
        else
            return x
        end
    end
end

function ExactDiagonalization.simplify(
    x::ParticleProjectorSumOperator{BR, S};
    tol::Real=Base.rtoldefault(real(S))
) where {BR, S}
    PPUO = ParticleProjectorUnitOperator{BR, S}
    PPSO = ParticleProjectorSumOperator{BR, S}
    terms::Vector{PPUO} = filter(x -> !iszero(x), simplify.(x.terms))
    isempty(terms) && return zero(PPUO)

    sort!(terms; lt=isless)
    new_terms = PPUO[]

    t1 = terms[1]

    bm::BR = t1.bitmask
    br::BR = t1.bitrow
    bc::BR = t1.bitcol
    pbm::BR= t1.parity_bitmask
    am::S  = t1.amplitude

    for t in terms[2:end]
        if (bm == t.bitmask) && (br == t.bitrow) && (bc == t.bitcol) && (pbm == t.parity_bitmask)
            am += t.amplitude
        else
            if !isapprox(am, zero(S); atol=tol)
                push!(new_terms, PPUO(bm, br, bc, pbm, am))
            end
            bm = t.bitmask
            br = t.bitrow
            bc = t.bitcol
            am = t.amplitude
        end
    end

    if !isapprox(am, zero(S); atol=tol)
        push!(new_terms, ParticleProjectorUnitOperator(bm, br, bc, am))
    end

    if isempty(new_terms)
        return NullOperator()
    end

    R = real(S)
    if S <: Complex && isapprox(maximum(abs(imag(t.amplitude)) for t in new_terms), zero(R); atol=tol)
        if length(new_terms) == 1
            return real(new_terms[1])
        else
            return ParticleProjectorSumOperator(real.(new_terms))
        end
    else
        if length(new_terms) == 1
            return new_terms[1]
        else
            return ParticleProjectorSumOperator(new_terms)
        end
    end
end
