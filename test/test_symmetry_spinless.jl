using Test
using LatticeTools
using QuantumHamiltonian
using QuantumHamiltonianParticle
using LinearAlgebra

@testset "Symmetry" begin
    f = Fermion("f")
    particle_sector = ParticleSector(f)
    site = ParticleSite([
        ParticleState(particle_sector, "_", [0], ( 0)),
        ParticleState(particle_sector, "f", [1], ( 1)),
    ])

    nsites = 4
    unitcell = makeunitcell(1.0; SiteType=String)
    addsite!(unitcell, "fermion orbital", FractCoord([0], [0.0]))
    lattice = makelattice(unitcell, nsites)

    tsym = FiniteTranslationSymmetry(lattice)
    psym = project(PointSymmetryDatabase.get3d(2), [1 0 0;])

    ssym = SymmorphicSymmetry(tsym, psym)
    ssymbed = embed(lattice, ssym)

    hilbert_space = ParticleHilbertSpace([site for i in 1:nsites])
    f_dag(i::Integer) = ParticleLadderUnit(particle_sector, 1, i, CREATION)
    f(i::Integer) = ParticleLadderUnit(particle_sector, 1, i, ANNIHILATION)

    hopping_hamiltonian = embed(hilbert_space,
        simplify(sum(
            let j = mod(i, nsites) + 1
                f_dag(i) * f(j) + f_dag(j) * f(i)
            end for i in 1:nsites
        ))
    )

    interaction_hamiltonian = embed(hilbert_space,
        simplify(sum(
            let j = mod(i, nsites) + 1
                f_dag(i) * f(i) * f_dag(j) * f(j) * 10
            end for i in 1:nsites
        ))
    )

    hamiltonian = make_projector_operator(hopping_hamiltonian + interaction_hamiltonian)
    hsr = represent(hilbert_space)
    rhsr = symmetry_reduce(hsr, first(get_irrep_components(ssymbed)))
    @testset "reduced hilbert space" begin
        @test get_basis_list(rhsr) == UInt[0b0000, 0b0001]
        # 0000 0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111
        @test rhsr.basis_mapping_index == [
            1,    2,   2,  -1,   2,  -1,  -1,  -1,   2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        ]
        @test rhsr.basis_mapping_amplitude[1] ≈ 1.0
        @test rhsr.basis_mapping_amplitude[2] ≈ 0.5
        @test rhsr.basis_mapping_amplitude[3] ≈ 0.5
        @test rhsr.basis_mapping_amplitude[5] ≈ 0.5
        @test rhsr.basis_mapping_amplitude[9] ≈ 0.5
    end

    @testset "eigenvalues" begin
        v1 = eigvals(Hermitian(Matrix(represent(hsr, hamiltonian))))
        v2 = Float64[]
        for sic in get_irrep_components(ssymbed)
            rhsr = symmetry_reduce(hsr, sic)
            h = represent(rhsr, hamiltonian)
            w = eigvals(Hermitian(Matrix(h)))
            append!(v2, w)
        end
        sort!(v2)
        @test maximum(abs2.(v1 - v2)) < 1E-8
    end
end
