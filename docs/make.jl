using Survival, Documenter, StatsBase, StatsAPI, AlgebraOfGraphics, CairoMakie, RDatasets

makedocs(
    modules = [Survival, Base.get_extension(Survival, :SurvivalAlgebraOfGraphicsExt)::Module],
    sitename = "Survival.jl",
    authors = "Alex Arslan",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Event Times" => "events.md",
        "Kaplan-Meier" => "km.md",
        "Nelson-Aalen" => "na.md",
        "Cox" => "cox.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaStats/Survival.jl.git",
    target = "build",
)
