using Particle

ElectronUp = Fermion{Symbol("e↑")}
ElectronDown = Fermion{Symbol("e↓")}

em = ParticleState("em", (ElectronUp(0), ElectronDown(0)), (0, 0))
up = ParticleState("up", (ElectronUp(1), ElectronDown(0)), (1, 1))
dn = ParticleState("dn", (ElectronUp(0), ElectronDown(1)), (1,-1))
ud = ParticleState("ud", (ElectronUp(1), ElectronDown(1)), (2, 0))

hubbard_states = [em, up, dn, ud]
tj_states = [em, up, dn]

hubbard_site = ParticleSite([em, up, dn, ud])
tj_site = ParticleSite([em, up, dn])
spin_site = ParticleSite([up, dn])

hs = ParticleHilbertSpace([tj_site, spin_site, tj_site])

for idx in CartesianIndices((1:3, 1:2, 1:3))
    s = join([site.states[i].name for (site, i) in zip(hs.sites, idx.I)], "-")
    l = getparticlelocationlist(hs, idx)
    println(idx.I, " | ", s, " | ", l)
end

ParticleSectorType = ParticleSector{Tuple{ElectronUp, ElectronDown}}
cupdag(i::Integer) = ParticleCreationOperator{ParticleSectorType}(1, i)
cdndag(i::Integer) = ParticleCreationOperator{ParticleSectorType}(2, i)
cup(i::Integer) = ParticleAnnihilationOperator{ParticleSectorType}(1, i)
cdn(i::Integer) = ParticleAnnihilationOperator{ParticleSectorType}(2, i)

@show shortstring( cupdag(3) * cup(3) )
@show shortstring( cupdag(3) * cdn(1) )

@show getstateindex(tj_site, ParticleSector(ElectronUp(0), ElectronDown(0)))
@show getstateindex(tj_site, ParticleSector(ElectronUp(1), ElectronDown(0)))
@show getstateindex(tj_site, ParticleSector(ElectronUp(0), ElectronDown(1)))

@show getparticlebitstring(hs, 0x0, UInt)

@show apply(hs, cupdag(3), 0x0)

@show apply(hs, cdndag(3), 0x0)



#=
ElectronUp = make_particle(Fermion)
ElectronDn = make_particle(Fermion)
em = ParticleState("em", (0, 0), ElectronUp => 0, ElectronDown => 0)
up = ParticleState("up", (0, 1), ElectronUp => 1, ElectronDown => 0)
dn = ParticleState("dn", (1, 0), ElectronUp => 0, ElectronDown => 1)
ud = ParticleState("ud", (1, 1), ElectronUp => 1, ElectronDown => 1)
s1 = ParticleSite([em, up, dn, ud])
s2 = ParticleSite([em, ud])
=#

;
