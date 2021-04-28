var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [QuantumHamiltonianParticle]","category":"page"},{"location":"api/#QuantumHamiltonianParticle.ParticleHilbertSpace","page":"API","title":"QuantumHamiltonianParticle.ParticleHilbertSpace","text":"ParticleHilbertSpace{PS, BR, QN} <: AbstractHilbertSpace{QN}\n\nParticle Hilbert space.\n\nExample\n\nParticleHilbertSpace([site1, site2, site3, site4])\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantumHamiltonianParticle.ParticleSite","page":"API","title":"QuantumHamiltonianParticle.ParticleSite","text":"ParticleSite{PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}\n\nParticle site.\n\nExample\n\nParticleSite([state1, state2, state3, state4])\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantumHamiltonianParticle.ParticleState","page":"API","title":"QuantumHamiltonianParticle.ParticleState","text":"ParticleState{PS, BR, QN}\n\nA state, represented by particle occupation.\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantumHamiltonianParticle.ParticleState-Union{Tuple{BR}, Tuple{PS}, Tuple{Type{PS}, AbstractString, AbstractVector{var\"#s24\"} where var\"#s24\"<:Integer}, Tuple{Type{PS}, AbstractString, AbstractVector{var\"#s25\"} where var\"#s25\"<:Integer, Type{BR}}} where {PS<:ParticleSector, BR<:Unsigned}","page":"API","title":"QuantumHamiltonianParticle.ParticleState","text":"ParticleState(::Type{PS}, name, occvec, ::Type{BR}=UInt)\n\nCreate a particle state with no quantum number.\n\nArguments\n\nPS: particle sector\nname: name of the state\noccvec::AbstractVector{<:Integer}: occupation Vector\nBR: binary type\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonianParticle.ParticleState-Union{Tuple{QN}, Tuple{BR}, Tuple{PS}, Tuple{Type{PS}, AbstractString, AbstractVector{var\"#s18\"} where var\"#s18\"<:Integer, QN}, Tuple{Type{PS}, AbstractString, AbstractVector{var\"#s19\"} where var\"#s19\"<:Integer, QN, Type{BR}}} where {PS<:ParticleSector, BR<:Unsigned, QN<:Tuple{Vararg{Integer, N} where N}}","page":"API","title":"QuantumHamiltonianParticle.ParticleState","text":"ParticleState(::Type{PS}, name, occvec, quantum_number, ::Type{BR}=UInt)\n\nCreate a particle state\n\nArguments\n\nPS: particle sector\nname: name of the state\noccvec::AbstractVector{<:Integer}: occupation Vector\nquantum_number: quantum number (tuple or integer)\nBR: binary type\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonianParticle.Spin","page":"API","title":"QuantumHamiltonianParticle.Spin","text":"Spin{Species, N}\n\nSpin type. S⁺ is annihilation, and S⁻ is creation.\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantumHamiltonian.bitoffset-Union{Tuple{QN}, Tuple{BR}, Tuple{PS}, Tuple{ParticleHilbertSpace{PS, BR, QN}, Integer, Integer}} where {PS, BR, QN}","page":"API","title":"QuantumHamiltonian.bitoffset","text":"bitoffset(phs, iptl, isite)\n\nGet the bit offset of the particle iptl at site isite.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonian.bitoffset-Union{Tuple{QN}, Tuple{BR}, Tuple{PS}, Tuple{ParticleHilbertSpace{PS, BR, QN}, Integer}} where {PS, BR, QN}","page":"API","title":"QuantumHamiltonian.bitoffset","text":"bitoffset(phs, isite)\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonian.bitwidth-Tuple{ParticleHilbertSpace}","page":"API","title":"QuantumHamiltonian.bitwidth","text":"bitwidth(phs)\n\nReturn number of bits needed to represent basis states of phs.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonian.compress-Union{Tuple{BR2}, Tuple{QN}, Tuple{BR}, Tuple{PS}, Tuple{ParticleHilbertSpace{PS, BR, QN}, CartesianIndex}, Tuple{ParticleHilbertSpace{PS, BR, QN}, CartesianIndex, Type{BR2}}} where {PS, BR, QN, BR2<:Unsigned}","page":"API","title":"QuantumHamiltonian.compress","text":"compress(hs, indexarray, [type])\n\nReturn the binary representation of the basis state represented by indexarray, optionally in type type.\n\nArguments\n\nhs::ParticleHilbertSpace{PS, BR, QN}\nindexarray::CartesianIndex\ntype::Type{BR2}=BR\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonian.extract-Tuple{ParticleHilbertSpace, Unsigned}","page":"API","title":"QuantumHamiltonian.extract","text":"extract(hs, occbin)\n\nReturn the CartesianIndex representation of the basis state represented by occbin.\n\nArguments\n\nhs::ParticleHilbertSpace\noccbin::Unsigned\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonian.get_bitmask-Union{Tuple{QN}, Tuple{BR}, Tuple{PS}, Tuple{ParticleHilbertSpace{PS, BR, QN}, Integer, Integer}} where {PS, BR, QN}","page":"API","title":"QuantumHamiltonian.get_bitmask","text":"get_bitmask(phs, [iptl, isite])\n\nGet the bit mask for the particles iptl at sites isite. iptl or isite can either be integer, a vector of integers, or colon :. Bitwise or is taken over list of iptl.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonian.get_quantum_number-Tuple{ParticleHilbertSpace, AbstractVector{var\"#s26\"} where var\"#s26\"<:Integer}","page":"API","title":"QuantumHamiltonian.get_quantum_number","text":"get_quantum_number(phs, statevec)\n\nGet quantum number of the basis state statevec.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonian.get_quantum_number-Tuple{ParticleHilbertSpace, Unsigned}","page":"API","title":"QuantumHamiltonian.get_quantum_number","text":"get_quantum_number(phs, bitrep)\n\nGet quantum number of the basis state bitrep.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonian.uncompress-Tuple{ParticleHilbertSpace, Unsigned}","page":"API","title":"QuantumHamiltonian.uncompress","text":"uncompress(hs, occbin)\n\nReturn the CartesianIndex representation of the basis state represented by occbin.\n\nArguments\n\nhs::ParticleHilbertSpace\noccbin::Unsigned\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonianParticle.electron_system-Tuple{}","page":"API","title":"QuantumHamiltonianParticle.electron_system","text":"electron_system()\n\nCreate an electron particle sector and creation/annihilation operators. Returns (particle_sector, c, cdag).candcdagare functions that take site index and spin, e.g.c(3, :up),cdag(4, :↓)`\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonianParticle.get_fermion_parity-Union{Tuple{PS}, Tuple{ParticleHilbertSpace, ParticleLadderUnit{PS, var\"#s27\", var\"#s26\"} where {var\"#s27\"<:Integer, var\"#s26\"<:Integer}, Unsigned}} where PS","page":"API","title":"QuantumHamiltonianParticle.get_fermion_parity","text":"get_fermion_parity(phs, op, bvec)\n\nGet the fermion parity (0 or 1) for when op is applied to bvec.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonianParticle.get_occupancy-Tuple{ParticleHilbertSpace, Integer, Integer, Unsigned}","page":"API","title":"QuantumHamiltonianParticle.get_occupancy","text":"get_occupancy(phs, iptl, isite, bvec::Unsigned)\n\nGet occupancy of particle iptl at site isite for the given basis state bvec.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonianParticle.get_parity_bitmask-Union{Tuple{QN}, Tuple{BR}, Tuple{PS}, Tuple{ParticleHilbertSpace{PS, BR, QN}, Integer, Integer}} where {PS, BR, QN}","page":"API","title":"QuantumHamiltonianParticle.get_parity_bitmask","text":"get_parity_bitmask(phs, iptl, isite)\n\nGet parity bitmask (i.e. position of the Wigner-Jordan string). Nonzero only for fermions.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonianParticle.locvec2statevec-Union{Tuple{QN}, Tuple{BR}, Tuple{PS}, Tuple{ParticleHilbertSpace{PS, BR, QN}, AbstractVector{var\"#s26\"} where var\"#s26\"<:(AbstractVector{var\"#s25\"} where var\"#s25\"<:Integer)}} where {PS, BR, QN}","page":"API","title":"QuantumHamiltonianParticle.locvec2statevec","text":"locvec2statevec\n\nparticles : Vector of (Vector of particle location)\n\nExample\n\n  P S 1 2 3 4 5 = |0, e↓, e↑, m, m⟩\n  e↑  0 0 1 0 0\n  e↓  0 1 0 0 0 = c†(m,4) c†(m,5) c†(e↑,3) c†(e↓,2) |Ω⟩\n  m   0 0 0 1 1\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumHamiltonianParticle.set_occupancy-Union{Tuple{BR2}, Tuple{QN}, Tuple{BR}, Tuple{PS}, Tuple{ParticleHilbertSpace{PS, BR, QN}, Integer, Integer, BR2, Integer}} where {PS, BR, QN, BR2}","page":"API","title":"QuantumHamiltonianParticle.set_occupancy","text":"set_occupancy(phs, iptl, isite, bvec::Unsigned, count)\n\nSet occupancy of particle iptl at site isite for the given basis state bvec to count.\n\n\n\n\n\n","category":"method"},{"location":"particle/#QuantumHamiltonianParticle-Sector","page":"Particle","title":"QuantumHamiltonianParticle Sector","text":"","category":"section"},{"location":"particle/#Particle-Types","page":"Particle","title":"Particle Types","text":"","category":"section"},{"location":"particle/","page":"Particle","title":"Particle","text":"DocTestSetup = quote\n    using QuantumHamiltonianParticle\nend","category":"page"},{"location":"particle/","page":"Particle","title":"Particle","text":"julia> f = Fermion(\"f\")\nFermion{:f}()\n\njulia> b = Boson(\"m\", 10)\nBoson{:m, 10}()\n\njulia> h = HardcoreBoson(\"h\")\nHardcoreBoson{:h}()\n\njulia> s = Spin(\"s=1/2\", 2)\nSpin{Symbol(\"s=1/2\"), 2}()\n\njulia> isboson(f), isfermion(f), isspin(f)\n(false, true, false)\n\njulia> isboson(b), isfermion(b), isspin(b)\n(true, false, false)\n\njulia> isboson(h), isfermion(h), isspin(h)\n(true, false, false)\n\njulia> isboson(s), isfermion(s), isspin(s)\n(false, false, true)\n\njulia> exchangesign(f), exchangesign(b), exchangesign(h), exchangesign(s)\n(-1, 1, 1, 1)\n\njulia> maxoccupancy(f), maxoccupancy(b), maxoccupancy(h), maxoccupancy(s)\n(1, 10, 1, 2)\n\njulia> using QuantumHamiltonian\n\njulia> bitwidth(f), bitwidth(b), bitwidth(h), bitwidth(s)\n(1, 4, 1, 2)","category":"page"},{"location":"particle/#Particle-Sector","page":"Particle","title":"Particle Sector","text":"","category":"section"},{"location":"particle/","page":"Particle","title":"Particle","text":"julia> electron_up, electron_dn = Fermion(\"e↑\"), Fermion(\"e↓\")\n(Fermion{Symbol(\"e↑\")}(), Fermion{Symbol(\"e↓\")}())\n\njulia> particle_sector = ParticleSector(electron_up, electron_dn)\nParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}()\n\njulia> numspecies(particle_sector)\n2\n\njulia> speciescount(particle_sector)\n2\n\njulia> getspecies(particle_sector)\n(Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")})\n\njulia> getspecies(particle_sector, 1)\nFermion{Symbol(\"e↑\")}\n\njulia> getspeciesname(particle_sector, 1)\nSymbol(\"e↑\")\n\njulia> exchangesign(particle_sector, 1)\n-1\n\njulia> using QuantumHamiltonian\n\njulia> bitwidth(particle_sector)\n2\n\njulia> bitwidth(particle_sector, 2)\n1\n\njulia> bitoffset(particle_sector, 2)\n1\n\njulia> bitoffset(particle_sector)\n3-element Vector{Int64}:\n 0\n 1\n 2","category":"page"},{"location":"particle/","page":"Particle","title":"Particle","text":"DocTestSetup = quote\n    using QuantumHamiltonianParticle\n    using QuantumHamiltonian\n    electron_up, electron_dn = Fermion(\"e↑\"), Fermion(\"e↓\")\n    particle_sector = ParticleSector(electron_up, electron_dn)\nend","category":"page"},{"location":"particle/","page":"Particle","title":"Particle","text":"julia> PS = typeof(particle_sector)\nParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}\n\njulia> get_bitmask(particle_sector, 2)\n0x0000000000000002\n\njulia> get_bitmask(PS, 2)\n0x0000000000000002\n\njulia> compress(particle_sector, [0,1])\n0x0000000000000002\n\njulia> compress(particle_sector, [1,1])\n0x0000000000000003\n\njulia> extract(particle_sector, 0x2)\n2-element Vector{Int64}:\n 0\n 1","category":"page"},{"location":"particle/","page":"Particle","title":"Particle","text":"DocTestSetup = nothing","category":"page"},{"location":"ladder/#Ladder","page":"Ladder","title":"Ladder","text":"","category":"section"},{"location":"ladder/","page":"Ladder","title":"Ladder","text":"DocTestSetup = quote\n    using QuantumHamiltonianParticle\n    electron_up, electron_dn = Fermion(\"e↑\"), Fermion(\"e↓\")\n    particle_sector = ParticleSector(electron_up, electron_dn)\n    cup(i) = ParticleLadderUnit(particle_sector, 1, i, ANNIHILATION)\n    cdn(i) = ParticleLadderUnit(particle_sector, 2, i, ANNIHILATION)\n    cupdag(i) = ParticleLadderUnit(particle_sector, 1, i, CREATION)\n    cdndag(i) = ParticleLadderUnit(particle_sector, 2, i, CREATION)\nend","category":"page"},{"location":"ladder/","page":"Ladder","title":"Ladder","text":"Set up","category":"page"},{"location":"ladder/","page":"Ladder","title":"Ladder","text":"using QuantumHamiltonianParticle\nelectron_up, electron_dn = Fermion(\"e↑\"), Fermion(\"e↓\")\nparticle_sector = ParticleSector(electron_up, electron_dn)\ncup(i) = ParticleLadderUnit(particle_sector, 1, i, ANNIHILATION)\ncdn(i) = ParticleLadderUnit(particle_sector, 2, i, ANNIHILATION)\ncupdag(i) = ParticleLadderUnit(particle_sector, 1, i, CREATION)\ncdndag(i) = ParticleLadderUnit(particle_sector, 2, i, CREATION)","category":"page"},{"location":"ladder/#Ladder-Unit","page":"Ladder","title":"Ladder Unit","text":"","category":"section"},{"location":"ladder/","page":"Ladder","title":"Ladder","text":"julia> cup(10)\nParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}, Int64, Int64}(1, 10, ANNIHILATION)\n\njulia> cup(10) == cup(10)\ntrue\n\njulia> cup(10) == cdn(10)\nfalse\n\njulia> exchangesign(cup(10), cup(3))\n-1\n\njulia> exchangesign(cup(10), cdn(3))\n1\n\njulia> maxoccupancy(cup(1))\n1\n\njulia> iszero(cup(1))\nfalse\n\njulia> adjoint(cup(10))\nParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}, Int64, Int64}(1, 10, CREATION)\n\njulia> using LinearAlgebra\n\njulia> ishermitian(cup(10))\nfalse\n\njulia> prettyprintln(cup(1))\nψ(e↑,1)","category":"page"},{"location":"ladder/#Ladder-Product","page":"Ladder","title":"Ladder Product","text":"","category":"section"},{"location":"ladder/","page":"Ladder","title":"Ladder","text":"DocTestSetup = quote\n    using QuantumHamiltonianParticle\n    electron_up, electron_dn = Fermion(\"e↑\"), Fermion(\"e↓\")\n    particle_sector = ParticleSector(electron_up, electron_dn)\n    cup(i) = ParticleLadderUnit(particle_sector, 1, i, ANNIHILATION)\n    cdn(i) = ParticleLadderUnit(particle_sector, 2, i, ANNIHILATION)\n    cupdag(i) = ParticleLadderUnit(particle_sector, 1, i, CREATION)\n    cdndag(i) = ParticleLadderUnit(particle_sector, 2, i, CREATION)\nend","category":"page"},{"location":"ladder/","page":"Ladder","title":"Ladder","text":"julia> cupdag(3) * cup(1)\nParticleLadderProduct{ParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}, Int64, Int64}(ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}, Int64, Int64}[ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}, Int64, Int64}(1, 3, CREATION), ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}, Int64, Int64}(1, 1, ANNIHILATION)])\n\njulia> adjoint(cupdag(3) * cup(1))\nParticleLadderProduct{ParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}, Int64, Int64}(ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}, Int64, Int64}[ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}, Int64, Int64}(1, 1, CREATION), ParticleLadderUnit{ParticleSector{Tuple{Fermion{Symbol(\"e↑\")}, Fermion{Symbol(\"e↓\")}}}, Int64, Int64}(1, 3, ANNIHILATION)])\n\njulia> using LinearAlgebra\n\njulia> ishermitian(cupdag(3) * cup(1))\nfalse\n\njulia> ishermitian(cupdag(1) * cup(1))\ntrue\n\njulia> prettyprintln(cupdag(3)*cdn(1))\nψ†(e↑,3)⋅ψ(e↓,1)","category":"page"},{"location":"#QuantumHamiltonianParticle","page":"Home","title":"QuantumHamiltonianParticle","text":"","category":"section"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"QuantumHamiltonianParticle.jl is an extension to QuantumHamiltonian.jl.","category":"page"}]
}
