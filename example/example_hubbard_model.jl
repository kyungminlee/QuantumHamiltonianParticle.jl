using LatticeTools
using ExactDiagonalization
using Particle

fermion = Fermion{:c}()
particle_sector = make_particle_sector(fermion)


c_dag(i, sp) = LadderUnitOperator(1, 2*(i-1) + sp, CREATION)
c(i, sp) = LadderUnitOperator(1, 2*(i-1) + sp, ANNIHILATION)

site = ParticleSite([
    ParticleState(particle_sector, "_", [0], (0,)),
    ParticleState(particle_sector, "f", [1], (1,)),
])

nsites = 4
hs = ParticleHilbertSpace([site for i in 1:nsites for sp in 1:2])

hopping = sum(
    let j = mod(i, nsites) + 1
        c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
    end
        for i in 1:nsites
        for sp in 1:2
)
hopping2 = make_projector_operator(hs, hopping)

interaction = sum(
    c_dag(i, 1) * c(i, 1) * c_dag(i, 2) * c(i, 2)
    for i in 1:nsites
)

t = 1.0
U = 1.0
hamiltonian = -t * hopping + U * interaction
hamiltonian_proj = make_projector_operator(hs, hamiltonian)
#@show make_projector_operator(hs, c_dag(3,1)  * c_dag(4,2) * c(4,2) * c(3, 1))
matrix = zeros(Float64, (2^(nsites*2), 2^(nsites*2)))
for bcol in UInt(0):UInt(1<<(2*nsites)-1)
    iter = get_column_iterator(hamiltonian_proj, bcol)
    for (brow, ampl) in iter
        matrix[brow+1, bcol+1] += ampl
    end
end
