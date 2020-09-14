export prettyprint
import ExactDiagonalization.prettyprintln

function prettyprint(arg::LadderUnitOperator)
    print("ψ")
    arg.ladder == CREATION && print("†")
    print("(", particle_species_name(arg.particle_index), ",", arg.orbital, ")")
end

function prettyprint(arg::LadderProductOperator)
    first = true
    for f in arg.factors
        if !first
            print("⋅")
        end
        prettyprint(f)
        first = false
    end
end

function prettyprint(arg::LadderSumOperator)
    if isempty(arg)
        print("0")
    else
        t, a = arg.terms[1]
        print("(", a, ")")
        if !isempty(t)
            print("⋅")
            prettyprint(t)
        end
        for (t, a) in arg.terms[2:end]
            print(" + (", a, ")")
            if !isempty(t)
                print("⋅")
                prettyprint(t)
            end
        end
    end
end


# function prettyprintln(io::IO, arg::ParticleProjectionUnitOperator{BR}, prefix::AbstractString="") where {BR}
#     pad = sizeof(BR)*8
#     println(io, prefix, "ParticleProjectionUnitOperator")
#     println(io, prefix, "  - bitmask: ", string(arg.bitmask, base=2, pad=pad))
#     println(io, prefix, "  - bitrow : ", string(arg.bitrow,  base=2, pad=pad))
#     println(io, prefix, "  - bitcol : ", string(arg.bitcol,  base=2, pad=pad))
#     println(io, prefix, "  - pmask  : ", string(arg.pmask,   base=2, pad=pad))
#     println(io, prefix, "  - prow   : ", string(arg.prow,    base=2, pad=pad))
#     println(io, prefix, "  - pcol   : ", string(arg.pcol,    base=2, pad=pad))
#     println(io, prefix, "  - pcheck : ", string(arg.pcheck,  base=2, pad=pad))
# end

# function prettyprintln(io::IO, arg::ParticleProjectionSumOperator{BR, S}, prefix::AbstractString="") where {BR, S}
#     println(io, prefix, "ParticleProjectionSumOperator")
#     for (term, ampl) in arg.terms
#         println(io, prefix, "amplit : ", ampl)
#         prettyprintln(io, term, prefix*"  ")
#     end
# end


function prettyprintln(io::IO, arg::ParticleProjectionUnitOperator{BR}, prefix::AbstractString="") where {BR}
    pad = sizeof(BR)*8
    println(io, prefix, "ParticleProjectionUnitOperator")
    println(io, prefix, "  - bitmask       : ", string(arg.bitmask, base=2, pad=pad))
    println(io, prefix, "  - bitrow        : ", string(arg.bitrow,  base=2, pad=pad))
    println(io, prefix, "  - bitcol        : ", string(arg.bitcol,  base=2, pad=pad))
    println(io, prefix, "  - paritybitmask : ", string(arg.parity_bitmask, base=2, pad=pad))
    println(io, prefix, "  - amplitude     : ", arg.amplitude)
end

function prettyprintln(io::IO, arg::ParticleProjectionSumOperator{BR, S}, prefix::AbstractString="") where {BR, S}
    println(io, prefix, "ParticleProjectionSumOperator")
    for term in arg.terms
        prettyprintln(io, term, prefix*"  ")
    end
end
