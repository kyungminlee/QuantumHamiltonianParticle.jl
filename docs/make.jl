using Documenter
using Particle

makedocs(
    modules=[Particle],
    doctest=true,
    sitename="Particle.jl",
    format=Documenter.HTML(prettyurls=!("local" in ARGS)),
    authors="Kyungmin Lee",
    checkdocs=:all,
    pages = [
        "Home" => "index.md",
        "Particle" => "particle.md",
        "Ladder" => "ladder.md",
        "API" => "api.md",
    ]
)

deploydocs(repo="github.com/kyungminlee/Particle.jl.git", devbranch = "dev")
