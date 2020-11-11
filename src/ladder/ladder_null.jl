export ParticleLadderNull
import LinearAlgebra

struct ParticleLadderNull{PS<:ParticleSector} <:AbstractParticleLadder{PS, Int}
    ParticleLadderNull{PS}() where {PS} = new{PS}()
    ParticleLadderNull(::Type{PS}) where {PS} = new{PS}()
    ParticleLadderNull(::PS) where {PS} = new{PS}()
end

exchangesign(::ParticleLadderNull{PS}, ::ParticleLadderNull{PS}) where {PS} = 1
exchangesign(::ParticleLadderNull{PS}, ::AbstractParticleLadder{PS, S}) where {PS, S} = 1
exchangesign(::AbstractParticleLadder{PS, S}, ::ParticleLadderNull{PS}) where {PS, S} = 1
Base.iszero(::ParticleLadderNull) = true
Base.adjoint(x::ParticleLadderNull) = x
LinearAlgebra.ishermitian(::ParticleLadderNull) = true
