using Documenter, SimilaritySolver

makedocs(
    sitename = "SimilaritySolver.jl",
    modules = [SimilaritySolver],
    warnonly = true,
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "API Reference" => "api.md",
        "Algorithm" => "algorithm.md",
        "Architecture" => "architecture.md",
        "Comparison" => "comparison.md",
        "Known Limitations" => "limitations.md",
        "Contributing" => "contributing.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://elvispy.github.io/SimilaritySolver.jl",
    ),
)

deploydocs(
    repo = "github.com/elvispy/SimilaritySolver.jl.git",
    devbranch = "main",
)
