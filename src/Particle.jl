module Particle

include("util.jl")

include("particle_type.jl")
include("particle_sector.jl")
include("particle_algebra.jl")

include("particle_site.jl")
include("particle_hilbert.jl")
include("particle_repconv.jl")

include("particle_simplify.jl")

# Free representation
include("particle_projection_operator.jl")
include("particle_sum_operator.jl")

include("prettyprint.jl")

end # module
