using LatticeTools
using QuantumHamiltonian
using QuantumHamiltonianParticle

using LinearAlgebra

fermion_up = Fermion{:up}()
fermion_dn = Fermion{:dn}()
particle_sector = ParticleSector(fermion_up, fermion_dn)

c_dag(i, sp) = ParticleLadderUnit(particle_sector, sp, i, CREATION)
c(i, sp) = ParticleLadderUnit(particle_sector, sp, i, ANNIHILATION)

site = ParticleSite([
    ParticleState(particle_sector, "em", [0,0], (0, 0)),
    ParticleState(particle_sector, "up", [1,0], (1, 1)),
    ParticleState(particle_sector, "dn", [0,1], (1,-1)),
    ParticleState(particle_sector, "ud", [1,1], (2, 0)),
])

nsites = 4
hs = ParticleHilbertSpace([site for i in 1:nsites])



hopping_nn = sum(
    let j = mod(i, nsites)+1
        c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
    end
        for i in 1:nsites
        for sp in 1:2
)
hopping_2nn = sum(
    let j = mod(i+1, nsites)+1
        c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
    end
        for i in 1:nsites
        for sp in 1:2
)

interaction = sum(
    c_dag(i, 1) * c(i, 1) * c_dag(i, 2) * c(i, 2)
    for i in 1:nsites
)


t = 1.0
tp = 0.6
U = 8.0
hamiltonian = embed(hs, -t * hopping_nn - tp * hopping_2nn + U * interaction)

hsr = represent(hs)
hamrep = represent(hsr, hamiltonian)


unitcell = makeunitcell(1.0; SiteType=String)
addsite!(unitcell, "A", FractCoord([0], [0.0]))

lattice = makelattice(unitcell, 4)

tsym = TranslationSymmetry(lattice)
tsymbed = embed(lattice, tsym)

for tic in get_irrep_components(tsymbed)
    for (symop, ampl) in get_irrep_iterator(tic)

    end
end

#tsym = TranslationSymmetry(4*ones(Int, (1,1)))
