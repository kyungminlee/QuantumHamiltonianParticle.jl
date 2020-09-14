module Particle

using ExactDiagonalization

include("util.jl")

include("basic/particle_type.jl")
include("basic/particle_sector.jl")

include("hilbert/particle_site.jl")
include("hilbert/particle_hilbert.jl")
include("hilbert/particle_repconv.jl")

#include("ladder/ladder_algebra.jl")
include("ladder/ladder_unit.jl")
include("ladder/ladder_product.jl")
include("ladder/ladder_sum.jl")
include("ladder/ladder_parity.jl")
include("ladder/ladder_apply.jl")

include("simplify/particle_simplify.jl")

# Free representation
#include("projection/particle_projection_unit_operator.jl")
#include("projection/particle_projection_sum_operator.jl")

include("projection/pure_operator.jl")
include("projection/sum_operator.jl")
include("projection/operator_iterator.jl")


include("io/prettyprint.jl")

end # module
