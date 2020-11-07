export ParticleLadderNull
import LinearAlgebra

struct ParticleLadderNull{PS<:ParticleSector} <:AbstractParticleLadder{PS, Int}
    ParticleLadderNull(::Type{PS}) where {PS} = new{PS}()
    ParticleLadderNull{PS}() where {PS} = new{PS}()
end

exchangesign(::ParticleLadderNull, ::AbstractParticleLadder) = 1
exchangesign(::AbstractParticleLadder, ::ParticleLadderNull) = 1
Base.iszero(::ParticleLadderNull) = true
Base.adjoint(x::ParticleLadderNull) = x
LinearAlgebra.ishermitian(::ParticleLadderNull) = true
