using Test
using LatticeTools
using QuantumHamiltonian
using QuantumHamiltonianParticle
using LinearAlgebra

@testset begin
    f = Fermion("f")
    particle_sector = ParticleSector(f)
    site = ParticleSite([
        ParticleState(particle_sector, "_", [0], ( 0)),
        ParticleState(particle_sector, "f", [1], ( 1)),
    ])

    nsites = 6
    unitcell = makeunitcell(1.0; SiteType=String)
    addsite!(unitcell, "fermion orbital", FractCoord([0], [0.0]))
    lattice = makelattice(unitcell, nsites)

    tsym = TranslationSymmetry(lattice)
    psym = project(PointSymmetryDatabase.get3d(2), [1 0 0;])

    ssym = tsym ⋊ psym
    tsymbed = embed(lattice, tsym)
    ssymbed = embed(lattice, ssym)

    hilbert_space = ParticleHilbertSpace([site for i in 1:nsites])

    f_dag(i::Integer) = ParticleLadderUnit(particle_sector, 1, i, CREATION)
    f(i::Integer) = ParticleLadderUnit(particle_sector, 1, i, ANNIHILATION)

    hopping_hamiltonian = embed(hilbert_space,
        simplify(sum(
            let j = mod(i, nsites) + 1
                f_dag(i) * f(j) + f_dag(j) * f(i)
            end
                for i in 1:nsites
        ))
    )

    interaction_hamiltonian = embed(hilbert_space,
        simplify(sum(
            let j = mod(i, nsites) + 1
                f_dag(i) * f(i) * f_dag(j) * f(j)
            end
                for i in 1:nsites
        ))
    )
    # @show hopping_hamiltonian
    # @show interaction_hamiltonian
    hamiltonian = hopping_hamiltonian + interaction_hamiltonian
    hamiltonian = make_projector_operator(hamiltonian)

    hsr = represent(hilbert_space)
    v1 = eigvals(Hermitian(Matrix(represent(hsr, hamiltonian))))
    v1[abs2.(v1) .< Base.rtoldefault(Float64)] .= 0
    @show v1

    # @show isinvariant(hamiltonian)
    # @show hamiltonian

    v2 = Float64[]
    for sic in get_irrep_components(ssymbed)
        rhsr = symmetry_reduce(hsr, sic)
        h = represent(rhsr, hamiltonian)

        w = eigvals(Hermitian(Matrix(h)))
        append!(v2, w)
    end
    sort!(v2)

    @show v2

    dv = v1 - v2
    dv[abs2.(dv) .< 1E-16] .= 0
    @show dv

    #=
    for sic in get_irrep_components(ssymbed)
        rhsr = symmetry_reduce(hsr, sic)
        println("SIC: $(sic.normal.irrep_index) $(sic.rest.irrep_index)/$(sic.rest.irrep_component)")
        for (j, (i, a)) in enumerate(zip(rhsr.basis_mapping_index, rhsr.basis_mapping_amplitude))
            print(string(hsr.basis_list[j], base=2, pad=4), "\t")
            if i <= 0
                println(-1)
            else
                println(string(rhsr.basis_list[i], base=2, pad=4), "\t", a)
            end
        end
        println()
    end
    =#
end
