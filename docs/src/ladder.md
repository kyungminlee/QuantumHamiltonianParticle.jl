# Ladder

```@meta
DocTestSetup = quote
    using Particle
    electron_up, electron_dn = Fermion("e↑"), Fermion("e↓")
    particle_sector = ParticleSector(electron_up, electron_dn)
    cup(i) = ParticleLadderUnit(particle_sector, 1, i, ANNIHILATION)
    cdn(i) = ParticleLadderUnit(particle_sector, 2, i, ANNIHILATION)
    cupdag(i) = ParticleLadderUnit(particle_sector, 1, i, CREATION)
    cdndag(i) = ParticleLadderUnit(particle_sector, 2, i, CREATION)
end
```

Set up
```
using Particle
electron_up, electron_dn = Fermion("e↑"), Fermion("e↓")
particle_sector = ParticleSector(electron_up, electron_dn)
cup(i) = ParticleLadderUnit(particle_sector, 1, i, ANNIHILATION)
cdn(i) = ParticleLadderUnit(particle_sector, 2, i, ANNIHILATION)
cupdag(i) = ParticleLadderUnit(particle_sector, 1, i, CREATION)
cdndag(i) = ParticleLadderUnit(particle_sector, 2, i, CREATION)
```

## Ladder Unit

```jldoctest
julia> cup(10)
ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol("e↑")},Fermion{Symbol("e↓")}}},Int64,Int64}(1, 10, ANNIHILATION)

julia> cup(10) == cup(10)
true

julia> cup(10) == cdn(10)
false

julia> exchangesign(cup(10), cup(3))
-1

julia> exchangesign(cup(10), cdn(3))
1

julia> maxoccupancy(cup(1))
1

julia> iszero(cup(1))
false

julia> adjoint(cup(10))
ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol("e↑")},Fermion{Symbol("e↓")}}},Int64,Int64}(1, 10, CREATION)

julia> using LinearAlgebra

julia> ishermitian(cup(10))
false

julia> prettyprintln(cup(1))
ψ(e↑,1)
```

## Ladder Product

```@meta
DocTestSetup = quote
    using Particle
    electron_up, electron_dn = Fermion("e↑"), Fermion("e↓")
    particle_sector = ParticleSector(electron_up, electron_dn)
    cup(i) = ParticleLadderUnit(particle_sector, 1, i, ANNIHILATION)
    cdn(i) = ParticleLadderUnit(particle_sector, 2, i, ANNIHILATION)
    cupdag(i) = ParticleLadderUnit(particle_sector, 1, i, CREATION)
    cdndag(i) = ParticleLadderUnit(particle_sector, 2, i, CREATION)
end
```

```jldoctest
julia> cupdag(3) * cup(1)
ParticleLadderProduct{ParticleSector{Tuple{Fermion{Symbol("e↑")},Fermion{Symbol("e↓")}}},Int64,Int64}(ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol("e↑")},Fermion{Symbol("e↓")}}},Int64,Int64}[ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol("e↑")},Fermion{Symbol("e↓")}}},Int64,Int64}(1, 3, CREATION), ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol("e↑")},Fermion{Symbol("e↓")}}},Int64,Int64}(1, 1, ANNIHILATION)])

julia> adjoint(cupdag(3) * cup(1))
ParticleLadderProduct{ParticleSector{Tuple{Fermion{Symbol("e↑")},Fermion{Symbol("e↓")}}},Int64,Int64}(ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol("e↑")},Fermion{Symbol("e↓")}}},Int64,Int64}[ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol("e↑")},Fermion{Symbol("e↓")}}},Int64,Int64}(1, 1, CREATION), ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol("e↑")},Fermion{Symbol("e↓")}}},Int64,Int64}(1, 3, ANNIHILATION)])

julia> using LinearAlgebra

julia> ishermitian(cupdag(3) * cup(1))
false

julia> ishermitian(cupdag(1) * cup(1))
true

julia> prettyprintln(cupdag(3)*cdn(1))
ψ†(e↑,3)⋅ψ(e↓,1)
```
