using Survival, Documenter

makedocs(
    modules = [Survival],
    clean = false,
    format = :html,
    sitename = "Survival.jl",
    authors = "Alex Arslan",
    pages = [
        "Home" => "index.md",
        "Event Times" => "events.md",
        "Kaplan-Meier" => "km.md",
        "Nelson-Aalen" => "na.md",
        "Cox" => "cox.md",
    ],
)

deploydocs(
    repo = "github.com/ararslan/Survival.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
