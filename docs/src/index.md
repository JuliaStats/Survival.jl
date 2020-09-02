```@meta
DocTestSetup = :(using Survival, StatsBase)
CurrentModule = Survival
```

# Survival.jl

This package provides types and methods for performing
[survival analysis](https://en.wikipedia.org/wiki/Survival_analysis) in Julia.

## Installation

The package is registered in Julia's General package registry, and so it can 
be installed using `Pkg.add("Survival"))` or via the REPL package manager mode `pkg>add Survival`.

## Contents

```@contents
Pages = [
    "events.md",
    "km.md",
    "na.md",
    "cox.md",
]
Depth = 1
```
