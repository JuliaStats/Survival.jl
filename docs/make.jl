using Survival, Documenter

makedocs(
    modules = [Survival],
    clean = false,
    format = :html,
    sitename = "Survival.jl",
    authors = "Alex Arslan",
    pages = [
        "Home" => "index.md",
        "Kaplan-Meier" => "km.md",
    ],
)

deploydocs(
    repo = "github.com/ararslan/Survival.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
