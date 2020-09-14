module Particle

using ExactDiagonalization

include("util.jl")

include("basic/particlesector.jl")

include("hilbertspace/particlesite.jl")
include("hilbertspace/particlehilbertspace.jl")

include("operators/ladderoperator.jl")
include("operators/productoperator.jl")
include("operators/sumoperator.jl")

# include("basic/particle_type.jl")
# include("basic/particle_sector.jl")

# include("ladder/ladder_algebra.jl")

# include("hilbert/particle_site.jl")
# include("hilbert/particle_hilbert.jl")
# include("hilbert/particle_repconv.jl")

# include("simplify/particle_simplify.jl")

# # Free representation
# include("projection/particle_projection_unit_operator.jl")
# include("projection/particle_projection_sum_operator.jl")

# include("io/prettyprint.jl")

end # module
