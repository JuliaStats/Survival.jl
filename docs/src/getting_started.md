# Getting Started

## Cox Proportional Hazards model
This tutorial shows how to fit a cox model and set up the `EventTimes` for right
censored events.  

### Dataset
We first read in a data set containing criminal recidivism data, which can be
found in the RcmdrPlugin.survival package in R, with the following link.
[RcmdrPlugin.survival](https://www.rdocumentation.org/packages/RcmdrPlugin.survival/versions/1.2-2/topics/Rossi)  

The `rossi` data set is originally from Rossi et al. (1980). 
The data pertain to 432 convicts who were released from Maryland state prisons
in the 1970s and who were followed up for one year after release. 
Half the released convicts were assigned at random to an experimental treatment
in which they were given financial aid; half did not receive aid.

### Fitting a cox model
To read in the dataset in the `DataFrame` format we use a few packages

```julia
julia> using Survival
julia> using DataFrames
julia> using StatsModels
julia> using CSV

julia> rossi = CSV.read(abspath(joinpath(pwd(), "..", "..", "test", "data", "rossi.csv")), DataFrame);
```

We now add the event times in the format (times, non-censoring indicator) where
no censoring (arrest) is indicated by a `1`

```julia
julia> rossi.event = EventTime.(rossi.week, rossi.arrest .== 1)
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
To fit the cox model, use the `coxph` function

```julia
julia> fitcox = coxph(@formula(event ~ fin + age + race + wexp + mar + paro + prio), rossi)
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

### Utilities

Several functions act on a `CoxModel` object.  

We can extract coefficients and the variance matrix of the fitted parameter
vector.  

The estimates can be extracted by

```julia
julia> coef(fitcox)
7-element Vector{Float64}:
 -0.3794158823466362
 -0.057429889713653676
  0.313920393830735
 -0.14980964863737226
 -0.43372380447995285
 -0.08486148372086805
  0.09152099594619753
```
and standard errors by `stderror(fitcox)`.  
The entire variance matrix is extracted by

```julia
julia> vcov(fitcox)
7×7 Matrix{Float64}:
  0.0366261    -0.000258838  …   0.00266155   -0.00020223
 -0.000258838   0.000483947      0.000473736  -1.93851e-5
 -0.00374361   -0.00030333      -0.00189336    0.000504948
  6.55182e-5   -0.00143737      -0.00162298    0.00170044
  0.0018773    -0.00101833      -0.004449     -0.000483896
  0.00266155    0.000473736  …   0.0383206     0.000842758
 -0.00020223   -1.93851e-5       0.000842758   0.000820643
```

Finally there is a few functions that might be useful.  
`loglikelihood(fitcox)`   
`nullloglikelihood(fitcox)`   
`cox_fit.model.score`   
`dof(fitcox)`   
`nobs(fitcox)`   
`coeftable(fitcox).cols`   




