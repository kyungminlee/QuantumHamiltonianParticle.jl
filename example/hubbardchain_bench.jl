using LatticeTools
using QuantumHamiltonian
using QuantumHamiltonianParticle
using LinearAlgebra
using BenchmarkTools
#=
fermion_up = Fermion{:up}()
fermion_dn = Fermion{:dn}()
particle_sector = ParticleSector(fermion_up, fermion_dn)

c_dag(i, sp) = ParticleLadderUnit(particle_sector, sp, i, CREATION)
c(i, sp) = ParticleLadderUnit(particle_sector, sp, i, ANNIHILATION)
=#

ps, c, c_dag = electron_system()
site = ParticleSite([
    ParticleState(ps, "em", [0,0], (0, 0)),
    ParticleState(ps, "up", [1,0], (1, 0)),
    ParticleState(ps, "dn", [0,1], (0, 1)),
    ParticleState(ps, "ud", [1,1], (1, 1)),
])

nsites = 6
hs = ParticleHilbertSpace([site for i in 1:nsites])

hopping_nn = sum(
    let j = mod(i, nsites)+1
        c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
    end
        for i in 1:nsites
        for sp in [:up, :dn]
)
hopping_2nn = sum(
    let j = mod(i+1, nsites)+1
        c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
    end
        for i in 1:nsites
        for sp in [:up, :dn]
)

interaction = sum(
    c_dag(i, :up) * c(i, :up) * c_dag(i, :dn) * c(i, :dn)
    for i in 1:nsites
)


t = 1.0
tp = 0.6
U = 8.0
hamiltonian = embed(hs, -t * hopping_nn - tp * hopping_2nn + U * interaction)
hamiltonian_p = make_projector_operator(hs,  -t * hopping_nn - tp * hopping_2nn + U * interaction)

hsr = represent_dict(hs)
hsr = represent_array(hs)

hamrep = represent(hsr, hamiltonian)
hamrep_p = represent(hsr, hamiltonian_p)

v = rand(Float64, size(hamrep, 1))
v2 = zeros(Float64, size(hamrep, 1))
@benchmark mul!(v2, hamrep, v)
@benchmark mul!(v2, hamrep_p, v)



unitcell = makeunitcell(1.0; SiteType=String)
addsite!(unitcell, "A", FractCoord([0], [0.0]))

lattice = makelattice(unitcell, nsites)

tsym = FiniteTranslationSymmetry(lattice)
psym = project(PointSymmetryDatabase.find3d("-1"), [1 0 0;])

tsymbed = embed(lattice, tsym)
psymbed = embed(lattice, psym)

tii = 1

tsa = collect(get_irrep_iterator(IrrepComponent(tsymbed, tii)))
psymbed_little = little_symmetry(tsymbed, tii, psymbed)
pii = 1
pic = 1
psa = collect(get_irrep_iterator(IrrepComponent(psymbed_little, pii, pic)))

ssa = make_product_irrep(psa, tsa)

hsr = represent_dict(hs);
hsr = represnet_array(hs);

rhsr = symmetry_reduce(hsr, ssa);
h1 = represent(rhsr, hamiltonian);
h2 = represent(rhsr, hamiltonian_p);

v = rand(Float64, size(h1, 1));
v2 = zeros(Float64, size(h1, 1));
@benchmark mul!(v2, h1, v)
@benchmark mul!(v2, h2, v)

Profile.init()
Profile.clear()
@profile for i in 1:1000
    mul!(v2, h2, v)
end

ProfileView.view()
