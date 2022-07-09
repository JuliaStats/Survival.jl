using Survival, Documenter, StatsBase

makedocs(
    modules = [Survival],
    sitename = "Survival.jl",
    authors = "Alex Arslan",
    pages = [
        "Home" => "index.md",
        "Event Times" => "events.md",
        "Kaplan-Meier" => "km.md",
        "Nelson-Aalen" => "na.md",
        "Cox" => "cox.md",
        "Getting Started" => "getting_started.md"
    ],
)

deploydocs(
    repo = "github.com/JuliaStats/Survival.jl.git",
    target = "build",
)
