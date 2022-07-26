# Getting Started

## Cox proportional hazards regression model

This tutorial shows how to fit a Cox model and set up the [`EventTime`](@ref)
values for right censored events.

### Dataset

In this example we'll use a table containing criminal
[recidivism](https://en.wikipedia.org/wiki/Recidivism) data from Rossi et al. (1980).
The data pertain to 432 convicts who were released from Maryland state prisons
in the 1970s and who were followed for one year after release.
The released convicts were randomly assigned to receive or not receive financial aid
with equal probability.
The outcome of interest is the time from release to rearrest.

The dataset is available as a CSV file in this package's `test/data/` directory.
To load the data as a `DataFrame`, we'll use the [CSV](https://github.com/JuliaData/CSV.jl)
and [DataFrames](https://github.com/JuliaData/DataFrames.jl) packages.

```julia-repl
julia> using Survival, StatsModels, CSV, DataFrames

julia> rossi = CSV.read(joinpath(pkgdir(Survival), "test", "data", "rossi.csv")), DataFrame);
```

### Fitting the model

We now construct the event times used as the model response.
The `EventTime` constructor accepts a time value and an indicator of whether the value
was observed (`true`) or right-censored (`false`).
The times in this data frame are in the column `week` and the arrest status in `arrest`.
An `arrest` value of 1 indicates an observed event (arrest) and 0 indicates censoring.

```julia-repl
julia> rossi.event = EventTime.(rossi.week, rossi.arrest .== 1);

julia> first(rossi, 10)
5×10 DataFrame
 Row │ arrest  week   fin    age    race   wexp   mar    paro   prio   event
     │ Int64   Int64  Int64  Int64  Int64  Int64  Int64  Int64  Int64  EventTim…
─────┼───────────────────────────────────────────────────────────────────────────
   1 │      1     20      0     27      1      0      0      1      3  20
   2 │      1     17      0     18      1      0      0      1      8  17
   3 │      1     25      0     19      0      1      0      1     13  25
   4 │      0     52      1     23      1      1      1      1      1  52+
   5 │      0     52      0     19      0      1      0      1      3  52+
```

To fit the Cox model, we can use `fit(CoxModel, ...)` or the shorthand `coxph(...)`.

```julia-repl
julia> model = coxph(@formula(event ~ fin + age + race + wexp + mar + paro + prio), rossi)
StatsModels.TableRegressionModel{CoxModel{Float64}, Matrix{Float64}}

event ~ fin + age + race + wexp + mar + paro + prio

Coefficients:
────────────────────────────────────────────────
        Estimate  Std.Error    z value  Pr(>|z|)
────────────────────────────────────────────────
fin   -0.379416   0.191379   -1.98253     0.0474
age   -0.0574299  0.0219988  -2.61059     0.0090
race   0.31392    0.307995    1.01924     0.3081
wexp  -0.14981    0.212226   -0.705898    0.4803
mar   -0.433724   0.38187    -1.13579     0.2560
paro  -0.0848615  0.195756   -0.433505    0.6646
prio   0.091521   0.0286469   3.1948      0.0014
────────────────────────────────────────────────
```

### Accessing values

Many of the common functions for accessing model parameters used in packages such as
[GLM](https://github.com/JuliaStats/GLM.jl) are extended for use with `CoxModel`s.

For example, the model coefficient estimates can be extracted with `coef`:

```julia-repl
julia> coef(model)
7-element Vector{Float64}:
 -0.3794158823466362
 -0.057429889713653676
  0.313920393830735
 -0.14980964863737226
 -0.43372380447995285
 -0.08486148372086805
  0.09152099594619753
```

Similarly, the standard errors of the estimates are accessible with `stderror(model)`,
the full variance-covariance matrix with `vcov(model)`, and the Wald confidence 
intervals with `confint(model)`. 

Other available functions include:
- `loglikelihood`, the log likelihood of the fitted model
- `nullloglikelihood`, the log likelihood of the null model
- `dof`, degrees of freedom
- `nobs`, the number of observations used to fit the model
- `coeftable`, a table of coefficient names, estimates, standard errors, z-values, and p-values
