using Documenter
using QuantumHamiltonianParticle

makedocs(
    modules=[QuantumHamiltonianParticle],
    doctest=true,
    sitename="QuantumHamiltonianParticle.jl",
    format=Documenter.HTML(prettyurls=!("local" in ARGS)),
    authors="Kyungmin Lee",
    checkdocs=:all,
    pages = [
        "Home" => "index.md",
        "Particle" => "particle.md",
        "Ladder" => "ladder.md",
        "Fermion Parity" => "fermionparity.md",
        "API" => "api.md",
    ]
)

deploydocs(;
    repo="github.com/kyungminlee/QuantumHamiltonianParticle.jl.git",
    devbranch="dev"
)
