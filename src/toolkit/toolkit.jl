using QuantumHamiltonian

export electron_system


"""
    electron_system()

Create an electron particle sector and creation/annihilation operators.
Returns `(particle_sector, c, cdag)``. `c` and `cdag` are functions that take
site index and spin, e.g. `c(3, :up)`, `cdag(4, :↓)`
"""
function electron_system()
    electron_up, electron_dn = Fermion("e↑"), Fermion("e↓")
    particle_sector = ParticleSector(electron_up, electron_dn)
    function c(i::Integer, spin::Symbol)
        if spin in [:up, :Up, :UP, :↑]
            return ParticleLadderUnit(particle_sector, 1, i, ANNIHILATION)
        elseif spin in [:dn, :Dn, :DN, :down, :Down, :DOWN, :↓]
            return ParticleLadderUnit(particle_sector, 2, i, ANNIHILATION)
        else
            throw(ArgumentError("unsupported spin $spin"))
        end
    end
    function cdag(i::Integer, spin::Symbol)
        if spin in [:up, :Up, :UP, :↑]
            return ParticleLadderUnit(particle_sector, 1, i, CREATION)
        elseif spin in [:dn, :Dn, :DN, :down, :Down, :DOWN, :↓]
            return ParticleLadderUnit(particle_sector, 2, i, CREATION)
        else
            throw(ArgumentError("unsupported spin $spin"))
        end
    end    
    return (particle_sector, c, cdag)
end
