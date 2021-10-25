# QuantumHamiltonianParticle Sector

## Particle Types

```@meta
DocTestSetup = quote
    using QuantumHamiltonianParticle
end
```

```jldoctest
julia> f = Fermion("f");

julia> b = Boson("m", 10);

julia> h = HardcoreBoson("h");

julia> s = Spin("s=1/2", 2);

julia> isboson(f), isfermion(f), isspin(f)
(false, true, false)

julia> isboson(b), isfermion(b), isspin(b)
(true, false, false)

julia> isboson(h), isfermion(h), isspin(h)
(true, false, false)

julia> isboson(s), isfermion(s), isspin(s)
(false, false, true)

julia> exchangesign(f), exchangesign(b), exchangesign(h), exchangesign(s)
(-1, 1, 1, 1)

julia> maxoccupancy(f), maxoccupancy(b), maxoccupancy(h), maxoccupancy(s)
(1, 10, 1, 2)

julia> using QuantumHamiltonian

julia> bitwidth(f), bitwidth(b), bitwidth(h), bitwidth(s)
(1, 4, 1, 2)
```

## Particle Sector

```jldoctest
julia> electron_up, electron_dn = Fermion("e↑"), Fermion("e↓");

julia> particle_sector = ParticleSector(electron_up, electron_dn);

julia> numspecies(particle_sector)
2

julia> speciescount(particle_sector)
2

julia> getspecies(particle_sector)
(Fermion{Symbol("e↑")}, Fermion{Symbol("e↓")})

julia> getspecies(particle_sector, 1)
Fermion{Symbol("e↑")}

julia> getspeciesname(particle_sector, 1)
Symbol("e↑")

julia> exchangesign(particle_sector, 1)
-1

julia> using QuantumHamiltonian

julia> bitwidth(particle_sector)
2

julia> bitwidth(particle_sector, 2)
1

julia> bitoffset(particle_sector, 2)
1

julia> bitoffset(particle_sector)
3-element Vector{Int64}:
 0
 1
 2
```

```@meta
DocTestSetup = quote
    using QuantumHamiltonianParticle
    using QuantumHamiltonian
    electron_up, electron_dn = Fermion("e↑"), Fermion("e↓")
    particle_sector = ParticleSector(electron_up, electron_dn)
end
```

```jldoctest
julia> get_bitmask(particle_sector, 2)
0x0000000000000002

julia> get_bitmask(PS, 2)
0x0000000000000002

julia> compress(particle_sector, [0,1])
0x0000000000000002

julia> compress(particle_sector, [1,1])
0x0000000000000003

julia> extract(particle_sector, 0x2)
2-element Vector{Int64}:
 0
 1
```

```@meta
DocTestSetup = nothing
```