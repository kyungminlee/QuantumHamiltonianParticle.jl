# Compare two particle representation vs one particle representation

using LatticeTools
using ExactDiagonalization
using Particle

using LinearAlgebra

# ====== Spin up and Spin down separately =====

function make_hubbard_twospecies(t::Real, tp::Real, U::Real, nsites::Integer)
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

    hs = ParticleHilbertSpace([site for i in 1:nsites])

    hopping = sum(
        let j = mod(i, nsites) + 1
            c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
        end
            for i in 1:nsites
            for sp in 1:2
    )
    hopping2 = sum(
        let j = mod(i+1, nsites) + 1
            c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
        end
            for i in 1:nsites
            for sp in 1:2
    )

    interaction = sum(
        c_dag(i, 1) * c(i, 1) * c_dag(i, 2) * c(i, 2)
        for i in 1:nsites
    )

    hamiltonian = -t * hopping - tp * hopping + U * interaction
    hamiltonian_proj = make_projector_operator(hs, hamiltonian)

    #@show make_projector_operator(hs, c_dag(3,1)  * c_dag(4,2) * c(4,2) * c(3, 1))

    matrix = zeros(Float64, (2^(nsites*2), 2^(nsites*2)))
    for bcol in UInt(0):UInt(1<<(2*nsites)-1)
        iter = get_column_iterator(hamiltonian_proj, bcol)
        for (brow, ampl) in iter
            matrix[brow+1, bcol+1] += ampl
        end
    end

    #hss = HilbertSpaceSector(hs, (nsites, mod(nsites, 2) ))
    hss = HilbertSpaceSector(hs, [(nsites, sp) for sp in -nsites:2:nsites])
    hssr = represent(hss)

    h_rep = represent(hssr, hamiltonian_proj) # DOES IT TAKE FERMION SIGN INTO ACCOUNT?
    return Matrix(h_rep)
end



# ====== Spin up and Spin down together =====

function make_hubbard_onespecies(t::Real, tp::Real, U::Real, nsites::Integer)
    fermion = Fermion{:e}()
    particle_sector = ParticleSector(fermion)

    c_dag(i, sp) = ParticleLadderUnit(particle_sector, 1, 2*(i-1) + sp, CREATION)
    c(i, sp) = ParticleLadderUnit(particle_sector, 1, 2*(i-1) + sp, ANNIHILATION)

    site = ParticleSite([
        ParticleState(particle_sector, "_", [0,], (0,)),
        ParticleState(particle_sector, "e", [1,], (1,)),
    ])

    hs = ParticleHilbertSpace([site for i in 1:nsites for sp in 1:2])

    hopping = sum(
        let j = mod(i, nsites) + 1
            c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
        end
            for i in 1:nsites
            for sp in 1:2
    )
    hopping2 =  sum(
        let j = mod(i+1, nsites) + 1
            c_dag(i, sp)*c(j, sp) + c_dag(j, sp)*c(i, sp)
        end
            for i in 1:nsites
            for sp in 1:2
    )

    interaction = sum(
        c_dag(i, 1) * c(i, 1) * c_dag(i, 2) * c(i, 2)
        for i in 1:nsites
    )

    hamiltonian = -t * hopping - tp * hopping + U * interaction
    hamiltonian_proj = make_projector_operator(hs, hamiltonian)

    #@show make_projector_operator(hs, c_dag(3,1)  * c_dag(4,2) * c(4,2) * c(3, 1))

    matrix = zeros(Float64, (2^(nsites*2), 2^(nsites*2)))
    for bcol in UInt(0):UInt(1<<(2*nsites)-1)
        iter = get_column_iterator(hamiltonian_proj, bcol)
        for (brow, ampl) in iter
            matrix[brow+1, bcol+1] += ampl
        end
    end

    #hss = HilbertSpaceSector(hs, (nsites, mod(nsites, 2) ))
    hss = HilbertSpaceSector(hs, [nsites]  )
    hssr = represent(hss)

    h_rep = represent(hssr, hamiltonian_proj) # DOES IT TAKE FERMION SIGN INTO ACCOUNT?
    return Matrix(h_rep)
end


t = 1.0
tp = 0.1
U = 4.0
h1 = make_hubbard_onespecies(t, tp, U, 4)
h2 = make_hubbard_onespecies(t, tp, U, 4)
@show h1 ≈ h2
@show eigvals(Hermitian(h1)) ≈ eigvals(Hermitian(h2))
