using LatticeTools
using QuantumHamiltonian
using QuantumHamiltonianParticle

using LinearAlgebra

boson_up = HardcoreBoson{:up}()
boson_dn = HardcoreBoson{:dn}()
particle_sector = ParticleSector(boson_up, boson_dn)

c_dag(i, sp) = ParticleLadderUnit(particle_sector, sp, i, CREATION)
c(i, sp) = ParticleLadderUnit(particle_sector, sp, i, ANNIHILATION)


site = ParticleSite([
    ParticleState(particle_sector, "em", [0,0], (0, 0)),
    ParticleState(particle_sector, "up", [1,0], (1, 1)),
    ParticleState(particle_sector, "dn", [0,1], (1,-1)),
    ParticleState(particle_sector, "ud", [1,1], (2, 0)),
])

nsites = 6
hs = ParticleHilbertSpace([site for i in 1:nsites])

hopping_nn = sum(
    let j = i+1
        c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
    end
        for i in 1:nsites-1
        for sp in 1:2
)

hopping_2nn = sum(
    let j = i+2
        c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
    end
        for i in 1:nsites-2
        for sp in 1:2
)

interaction = sum(
    c_dag(i, 1) * c(i, 1) * c_dag(i, 2) * c(i, 2)
    for i in 1:nsites
)

t = 1.0
tp = 0.6
U = 8.0
hamiltonian = -t * hopping_nn - tp * hopping_2nn + U * interaction
hamiltonian_proj = make_projector_operator(hs, hamiltonian)

matrix = zeros(Float64, (2^(nsites*2), 2^(nsites*2)))
for bcol in UInt(0):UInt(1<<(2*nsites)-1)
    iter = get_column_iterator(hamiltonian_proj, bcol)
    for (brow, ampl) in iter
        matrix[brow+1, bcol+1] += ampl
    end
end

println("HardcoreBoson")
@show minimum(eigvals(Hermitian(matrix)));
