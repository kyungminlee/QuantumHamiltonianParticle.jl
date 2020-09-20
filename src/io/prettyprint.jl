export prettyprint, prettyprintln
import ExactDiagonalization.prettyprintln

prettyprintln(xs...) = prettyprintln(stdout::IO, xs...)
prettyprint(xs...) = prettyprint(stdout::IO, xs...)

function prettyprint(io::IO, arg::LadderUnitOperator{PS, PI, OI}) where {PS, PI, OI}
    print(io, "ψ")
    arg.ladder == CREATION && print(io, "†")
    print(io, "(", getspeciesname(PS, arg.particle_index), ",", arg.orbital, ")")
end

function prettyprint(io::IO, arg::LadderProductOperator)
    first = true
    for f in arg.factors
        if !first
            print(io, "⋅")
        end
        prettyprint(io, f)
        first = false
    end
end

function prettyprint(io::IO, arg::LadderSumOperator)
    if isempty(arg.terms)
        print(io, "0")
    else
        t, a = arg.terms[1]
        print("(", a, ")")
        if !isempty(t.factors)
            print(io, "⋅")
            prettyprint(io, t)
        end
        for (t, a) in arg.terms[2:end]
            print(io, " + (", a, ")")
            if !isempty(t.factors)
                print(io, "⋅")
                prettyprint(io, t)
            end
        end
    end
end


function prettyprintln(io::IO, arg::LadderUnitOperator)
    prettyprint(io, arg)
    println(io)
end

function prettyprintln(io::IO, arg::LadderProductOperator)
    prettyprint(io, arg)
    println(io)
end

function prettyprintln(io::IO, arg::LadderSumOperator)
    prettyprint(io, arg)
    println(io)
end

# function prettyprintln(io::IO, arg::ParticleProjectorUnitOperator{BR}, prefix::AbstractString="") where {BR}
#     pad = sizeof(BR)*8
#     println(io, prefix, "ParticleProjectorUnitOperator")
#     println(io, prefix, "  - bitmask: ", string(arg.bitmask, base=2, pad=pad))
#     println(io, prefix, "  - bitrow : ", string(arg.bitrow,  base=2, pad=pad))
#     println(io, prefix, "  - bitcol : ", string(arg.bitcol,  base=2, pad=pad))
#     println(io, prefix, "  - pmask  : ", string(arg.pmask,   base=2, pad=pad))
#     println(io, prefix, "  - prow   : ", string(arg.prow,    base=2, pad=pad))
#     println(io, prefix, "  - pcol   : ", string(arg.pcol,    base=2, pad=pad))
#     println(io, prefix, "  - pcheck : ", string(arg.pcheck,  base=2, pad=pad))
# end

# function prettyprintln(io::IO, arg::ParticleProjectorSumOperator{BR, S}, prefix::AbstractString="") where {BR, S}
#     println(io, prefix, "ParticleProjectorSumOperator")
#     for (term, ampl) in arg.terms
#         println(io, prefix, "amplit : ", ampl)
#         prettyprintln(io, term, prefix*"  ")
#     end
# end


function prettyprintln(io::IO, arg::ParticleProjectorUnitOperator{BR}, prefix::AbstractString="") where {BR}
    pad = sizeof(BR)*8
    println(io, prefix, "ParticleProjectorUnitOperator")
    println(io, prefix, "  - bitmask       : ", string(arg.bitmask, base=2, pad=pad))
    println(io, prefix, "  - bitrow        : ", string(arg.bitrow,  base=2, pad=pad))
    println(io, prefix, "  - bitcol        : ", string(arg.bitcol,  base=2, pad=pad))
    println(io, prefix, "  - paritybitmask : ", string(arg.parity_bitmask, base=2, pad=pad))
    println(io, prefix, "  - amplitude     : ", arg.amplitude)
end

function prettyprintln(io::IO, arg::ParticleProjectorSumOperator{BR, S}, prefix::AbstractString="") where {BR, S}
    println(io, prefix, "ParticleProjectorSumOperator")
    for term in arg.terms
        prettyprintln(io, term, prefix*"  ")
    end
end
