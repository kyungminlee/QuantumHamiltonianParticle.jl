import ExactDiagonalization.apply
import ExactDiagonalization.get_row_iterator
import ExactDiagonalization.get_column_iterator


function get_row_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderUnitOperator{PS, <:Integer, <:Integer},
    bvec::BR,
    ::Type{S}=Float64,
) where {PS, BR, QN, S<:Number}
    particle = getspecies(PS, op.particle_index)
    occupancy_at_site = get_occupancy(hs, op.particle_index, op.orbital, bvec)

    newbvec::BR = zero(BR)
    match::Bool = false
    ampl::S = zero(S)

    if op.ladder == ANNIHILATION
        if occupancy_at_site < maxoccupancy(particle)
            if isfermion(particle)
                wj_parity = get_fermion_parity(hs, op, bvec)
                ampl = wj_parity == 0 ? one(S) : -one(S)
            elseif isboson(particle)
                ampl = Base.sqrt(S(occupancy_at_site+1))
            elseif isspin(particle)
                ampl = one(S)
            else
                throw(ArgumentError("unsupported particle type $particle"))
            end
            newbvec = set_occupancy(hs, op.particle_index, op.orbital, bvec, occupancy_at_site+1)
            match = true
        end
    else # op.ladder == CREATION
        if occupancy_at_site > 0
            if isfermion(particle)
                wj_parity = get_fermion_parity(hs, op, bvec)
                ampl = wj_parity == 0 ? one(S) : -one(S)
            elseif isboson(particle)
                ampl = Base.sqrt(S(occupancy_at_site))
            elseif isspin(particle)
                ampl = one(S)
            else
                throw(ArgumentError("unsupported particle type $particle"))
            end
            newbvec = set_occupancy(hs, op.particle_index, op.orbital, bvec, occupancy_at_site-1)
            match = true
        end
    end
    element::Pair{BR, S} = newbvec => ampl
    return (element for i in 1:(match ? 1 : 0))
end


function get_column_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderUnitOperator{PS, <:Integer, <:Integer},
    bvec::BR,
    ::Type{S}=Float64,
) where {PS, BR, QN, S<:Number}
    particle = getspecies(PS, op.particle_index)
    occupancy_at_site = get_occupancy(hs, op.particle_index, op.orbital, bvec)

    newbvec::BR = zero(BR)
    match::Bool = false
    ampl::S = zero(S)

    if op.ladder == CREATION
        if occupancy_at_site < maxoccupancy(particle)
            if isfermion(particle)
                wj_parity = get_fermion_parity(hs, op, bvec)
                ampl = wj_parity == 0 ? one(S) : -one(S)
            elseif isboson(particle)
                ampl = Base.sqrt(S(occupancy_at_site+1))
            elseif isspin(particle)
                ampl = one(S)
            else
                throw(ArgumentError("unsupported particle type $particle"))
            end
            newbvec = set_occupancy(hs, op.particle_index, op.orbital, bvec, occupancy_at_site+1)
            match = true
        end
    else # op.ladder == ANNIHILATION
        if occupancy_at_site > 0
            if isfermion(particle)
                wj_parity = get_fermion_parity(hs, op, bvec)
                ampl = wj_parity == 0 ? one(S) : -one(S)
            elseif isboson(particle)
                ampl = Base.sqrt(S(occupancy_at_site))
            elseif isspin(particle)
                ampl = one(S)
            else
                throw(ArgumentError("unsupported particle type $particle"))
            end
            newbvec = set_occupancy(hs, op.particle_index, op.orbital, bvec, occupancy_at_site-1)
            match = true
        end
    end
    element::Pair{BR, S} = newbvec => ampl
    return (element for i in 1:(match ? 1 : 0))
end


function get_row_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderProductOperator{PS, <:Integer, <:Integer},
    bvec::BR,
   ::Type{S}=Float64,
) where {PS, BR, QN, S<:Number}
    match::Bool = true
    ampl::S = one(S)
    for f in op.factors
        newout = get_row_iterator(hs, f, bvec, S)
        if isempty(newout) || iszero(ampl)
            bvec = zero(BR)
            ampl = zero(S)
            match = false
            break
        end
        (newbvec, newampl) = first(newout)
        @assert length(newout) == 1
        (newbvec, newampl) = first(newout)
        bvec = newbvec
        ampl *= newampl
    end
    element::Pair{BR, S} = (bvec => ampl)
    return (element for i in 1:(match ? 1 : 0))
end


function get_column_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderProductOperator{PS, <:Integer, <:Integer},
    bvec::BR,
    ::Type{S}=Float64,
) where {PS, BR, QN, S<:Number}
    match::Bool = true
    ampl::S = one(S)
    for f in reverse(op.factors)
        newout = get_column_iterator(hs, f, bvec, S)
        if isempty(newout) || iszero(ampl)
            bvec = zero(BR)
            ampl = zero(S)
            match = false
            break
        end
        (newbvec, newampl) = first(newout)
        @assert length(newout) == 1
        (newbvec, newampl) = first(newout)
        bvec = newbvec
        ampl *= newampl
    end
    element::Pair{BR, S} = (bvec => ampl)
    return (element for i in 1:(match ? 1 : 0))
end


function get_row_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderSumOperator{PS, <:Integer, <:Integer, S},
    bvec::BR,
    ::Type{Sout}=promote_type(S, Float64),
) where {PS, BR, QN, S, Sout}
    return (
        newbvec => a * newampl
            for (t, a) in op.terms
            for (newbvec, newampl) in get_row_iterator(hs, t, bvec, Sout)
    )
end


function get_column_iterator(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderSumOperator{PS, <:Integer, <:Integer, S},
    bvec::BR,
    ::Type{Sout}=promote_type(S, Float64),
) where {PS, BR, QN, S, Sout}
    return (
        newbvec => a * newampl
            for (t, a) in op.terms
            for (newbvec, newampl) in get_column_iterator(hs, t, bvec, Sout)
    )
end



#=
function apply(
    hs::ParticleHilbertSpace{PS, BR, QN},
    op::LadderUnitOperator{<:Integer, <:Integer},
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
    op::LadderProductOperator,
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
    op::LadderSumOperator{<:Integer, <:Integer, S},
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
