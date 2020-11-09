import ExactDiagonalization.apply
import ExactDiagonalization.get_row_iterator
import ExactDiagonalization.get_column_iterator
import ExactDiagonalization.get_element

function get_row_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderNull{PS},
    bvec::BR,
    ::Type{S}=Float64,
) where {PS, BR, QN, S}
    return ((zero(BR) => zero(S)) for i in 1:0)
end

function get_column_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderNull{PS},
    bvec::BR,
    ::Type{S}=Float64,
) where {PS, BR, QN, S}
    return ((zero(BR) => zero(S)) for i in 1:0)
end

function get_row_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderUnit{PS, <:Integer, <:Integer},
    bvec::BR,
    ::Type{S}=Float64,
) where {PS, BR, QN, S<:Number}
    particle = getspecies(PS, op.particle_index)
    occupancy_at_site = get_occupancy(hs, op.particle_index, op.orbital, bvec)
    new_occupancy_at_site = zero(occupancy_at_site)
    ref_occupancy_at_site = zero(occupancy_at_site)

    newbvec::BR = zero(BR)
    match::Bool = false
    ampl::S = zero(S)

    if op.ladder == ANNIHILATION
        if 0 <= occupancy_at_site < maxoccupancy(particle)
            new_occupancy_at_site = occupancy_at_site + 1
            ref_occupancy_at_site = occupancy_at_site + 1  # for amplitude sqrt(n) or sqrt(n+1)
            match = true
        end
    else # op.ladder == CREATION
        if 0 < occupancy_at_site <= maxoccupancy(particle)
            new_occupancy_at_site = occupancy_at_site - 1
            ref_occupancy_at_site = occupancy_at_site      # for amplitude sqrt(n) or sqrt(n+1)
            match = true
        end
    end
    if match
        if isfermion(particle)
            wj_parity = get_fermion_parity(hs, op, bvec)
            ampl = wj_parity == 0 ? one(S) : -one(S)
        elseif isboson(particle)
            ampl = Base.sqrt(S(ref_occupancy_at_site))
        elseif isspin(particle)
            ampl = one(S)
        else
            throw(ArgumentError("unsupported particle type $particle")) # COV_EXCL_LINE
        end
        newbvec = set_occupancy(hs, op.particle_index, op.orbital, bvec, new_occupancy_at_site)
    end
    element::Pair{BR, S} = newbvec => ampl
    return (element for i in 1:(match ? 1 : 0))
end


function get_column_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderUnit{PS, <:Integer, <:Integer},
    bvec::BR,
    ::Type{S}=Float64,
) where {PS, BR, QN, S<:Number}
    particle = getspecies(PS, op.particle_index)
    occupancy_at_site = get_occupancy(hs, op.particle_index, op.orbital, bvec)
    new_occupancy_at_site = zero(occupancy_at_site)
    ref_occupancy_at_site = zero(occupancy_at_site)

    newbvec::BR = zero(BR)
    match::Bool = false
    ampl::S = zero(S)

    if op.ladder == CREATION
        if 0 <= occupancy_at_site < maxoccupancy(particle)
            new_occupancy_at_site = occupancy_at_site + 1
            ref_occupancy_at_site = occupancy_at_site + 1   # for amplitude sqrt(n) or sqrt(n+1)
            match = true
        end
    else # op.ladder == ANNIHILATION
        if 0 < occupancy_at_site <= maxoccupancy(particle)
            new_occupancy_at_site = occupancy_at_site - 1
            ref_occupancy_at_site = occupancy_at_site       # for amplitude sqrt(n) or sqrt(n+1)
            match = true
        end
    end
    if match
        if isfermion(particle)
            wj_parity = get_fermion_parity(hs, op, bvec)
            ampl = wj_parity == 0 ? one(S) : -one(S)
        elseif isboson(particle)
            ampl = Base.sqrt(S(ref_occupancy_at_site))
        elseif isspin(particle)
            ampl = one(S)
        else
            throw(ArgumentError("unsupported particle type $particle")) # COV_EXCL_LINE
        end
        newbvec = set_occupancy(hs, op.particle_index, op.orbital, bvec, new_occupancy_at_site)
    end
    element::Pair{BR, S} = newbvec => ampl
    return (element for i in 1:(match ? 1 : 0))
end


function get_element(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderUnit{PS, <:Integer, <:Integer},
    brow::BR, bcol::BR,
    ::Type{S}=Float64,
) where {PS, BR, QN, S<:Number}
    particle = getspecies(PS, op.particle_index)
    occupancy_at_site = get_occupancy(hs, op.particle_index, op.orbital, bcol)
    new_occupancy_at_site = zero(occupancy_at_site)
    ref_occupancy_at_site = zero(occupancy_at_site)

    if op.ladder == CREATION
        if 0 <= occupancy_at_site < maxoccupancy(particle)
            new_occupancy_at_site = occupancy_at_site + 1
            ref_occupancy_at_site = occupancy_at_site + 1   # for amplitude sqrt(n) or sqrt(n+1)
        else
            return zero(S)
        end
    else
        if 0 < occupancy_at_site <= maxoccupancy(particle)
            new_occupancy_at_site = occupancy_at_site - 1
            ref_occupancy_at_site = occupancy_at_site       # for amplitude sqrt(n) or sqrt(n+1)
        else
            return zero(S)
        end
    end

    if brow == set_occupancy(hs, op.particle_index, op.orbital, bcol, new_occupancy_at_site)
        if isfermion(particle)
            wj_parity = get_fermion_parity(hs, op, bcol)
            ampl = wj_parity == 0 ? one(S) : -one(S)
        elseif isboson(particle)
            ampl = Base.sqrt(S(ref_occupancy_at_site))
        elseif isspin(particle)
            ampl = one(S)
        else
            throw(ArgumentError("unsupported particle type $particle")) # COV_EXCL_LINE
        end
        return ampl
    else
        return zero(S)
    end
end


function get_row_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderProduct{PS, <:Integer, <:Integer},
    bvec::BR,
   ::Type{S}=Float64,
) where {PS, BR, QN, S<:Number}
    match::Bool = true
    ampl::S = one(S)
    for f in op.factors
        newout = get_row_iterator(hs, f, bvec, S)
        if isempty(newout)
            bvec = zero(BR)
            ampl = zero(S)
            match = false
            break
        end
        @assert length(newout) == 1
        (newbvec, newampl) = first(newout)
        @assert !iszero(newampl)
        bvec = newbvec
        ampl *= newampl
    end
    element::Pair{BR, S} = (bvec => ampl)
    return (element for i in 1:(match ? 1 : 0))
end


function get_column_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderProduct{PS, <:Integer, <:Integer},
    bvec::BR,
    ::Type{S}=Float64,
) where {PS, BR, QN, S<:Number}
    match::Bool = true
    ampl::S = one(S)
    for f in reverse(op.factors)
        newout = get_column_iterator(hs, f, bvec, S)
        if isempty(newout)
            bvec = zero(BR)
            ampl = zero(S)
            match = false
            break
        end
        @assert length(newout) == 1
        (newbvec, newampl) = first(newout)
        @assert !iszero(newampl)
        bvec = newbvec
        ampl *= newampl
    end
    element::Pair{BR, S} = (bvec => ampl)
    return (element for i in 1:(match ? 1 : 0))
end


function get_element(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderProduct{PS, <:Integer, <:Integer},
    brow::BR, bcol::BR,
    ::Type{S}=Float64,
) where {PS, BR, QN, S<:Number}
    bvec = brow
    ampl::S = one(S)
    for f in op.factors
        newout = get_row_iterator(hs, f, bvec, S)
        isempty(newout) && return zero(S)
        @assert length(newout) == 1
        (newbvec, newampl) = first(newout)
        @assert !iszero(newampl)
        bvec = newbvec
        ampl *= newampl
    end
    return bvec == bcol ? ampl : zero(S)
end


function get_row_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderSum{PS, <:Integer, <:Integer, S},
    bvec::BR,
    ::Type{Sout}=promote_type(S, Float64),
) where {PS, BR, QN, S, Sout<:Number}
    return (
        newbvec => Sout(a * newampl)
            for (t, a) in op.terms
            for (newbvec, newampl) in get_row_iterator(hs, t, bvec, Sout)
    )
end


function get_column_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderSum{PS, <:Integer, <:Integer, S},
    bvec::BR,
    ::Type{Sout}=promote_type(S, Float64),
) where {PS, BR, QN, S, Sout<:Number}
    return (
        newbvec => Sout(a * newampl)
            for (t, a) in op.terms
            for (newbvec, newampl) in get_column_iterator(hs, t, bvec, Sout)
    )
end


function get_element(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderSum{PS, <:Integer, <:Integer, S},
    brow::BR, bcol::BR,
    ::Type{Sout}=promote_type(S, Float64),
) where {PS, BR, QN, S<:Number, Sout<:Number}
    return sum(get_element(hs, term, brow, bcol)*ampl for (term, ampl) in op.terms)::Sout
end



#=
function apply(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderUnit{<:Integer, <:Integer},
    bvec::BR2,
) where {PS, BR, QN, BR2<:Unsigned}

    particle = getspecies(PS, op.particle_index)
    occupancy_at_site = get_occupancy(hs, op.particle_index, op.orbital, bvec)
    # let
    #     occmat = occbin2occmat(hs, bvec)
    #     occupancy_at_site2 = occmat[op.particle_index, op.orbital]
    #     @assert occupancy_at_site == occupancy_at_site2
    # end
    if op.ladder == CREATION
        occupancy_at_site >= maxoccupancy(particle) && return [] #(~zero(bvec), 0.0)

        if isfermion(particle)
            wj_parity = get_fermion_parity(hs, op, bvec)
            # wj_parity2 = sum(occmat[op.particle_index, i] for i in 1:(op.orbital-1)) % 2
            # @show wj_parity, wj_parity2
            ampl = wj_parity == 0 ? 1.0 : -1.0
        elseif isboson(particle)
            ampl = sqrt(occupancy_at_site+1)
        else
            throw(ArgumentError("unsupported particle type $particle"))
        end
        # occmat[op.particle_index, op.orbital] += 1
        # newbvec = occmat2occbin(hs, occmat)
        newbvec = set_occupancy(hs, op.particle_index, op.orbital, bvec, occupancy_at_site+1)
        return [newbvec=>ampl for i in 1:1]
    else # op.ladder == ANNIHILATION
        occupancy_at_site <= 0 && return [] # (~zero(bvec), 0.0)

        if isfermion(particle)
            wj_parity = get_fermion_parity(hs, op, bvec)
            ampl = wj_parity == 0 ? 1.0 : -1.0
            # wj_parity2 = sum(occmat[op.particle_index, i] for i in 1:(op.orbital-1)) % 2
            # @show wj_parity, wj_parity2
        elseif isboson(particle)
            ampl = sqrt(occupation_at_site)
        else
            throw(ArgumentError("unsupported particle type $particle"))
        end
        # occmat[op.particle_index, op.orbital] -= 1
        # newbvec = occmat2occbin(hs, occmat)
        newbvec = set_occupancy(hs, op.particle_index, op.orbital, bvec, occupancy_at_site-1)
        return [newbvec=>ampl]
    end
end
=#

#=
function apply(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderProduct,
    bvec::BR2,
)::Vector{Pair{BR, Float64}} where {PS, BR, QN, BR2<:Unsigned}

    ampl = 1.0
    out = Pair{BR, Float64}[]
    for f in reverse(op.factors)
        newout = apply(hs, f, bvec)
        isempty(newout) && return Pair{BR, Float64}[]
        @assert length(newout) == 1
        (newbvec, newampl) = newout[1]
        ampl *= newampl
        bvec = newbvec
    end
    return [bvec=>ampl]
end
=#

#=
function apply(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::ParticleLadderSum{<:Integer, <:Integer, S},
    bvec::BR2,
)::Vector{Pair{BR, Float64}} where {PS, BR, QN, S, BR2<:Unsigned}
    Sout = promote_type(S, Float64)
    out = Pair{BR, Sout}[]
    return Pair{BR, Sout}[
        newbvec => a * newampl
            for (t, a) in reverse(op.terms)
            for (newbvec, newampl) in apply(hs, t, bvec)
    ]
    return out
end
=#
