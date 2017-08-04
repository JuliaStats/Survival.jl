import Survival
using ClobberingReload
using DataFrames
creload("Survival")
filepath = joinpath(Pkg.dir("Survival", "test"), "rossi.csv")
rossi = readtable(filepath)
rossi[:event] = Survival.EventTime.(rossi[:week],rossi[:arrest] .== 1)
Survival.coxph(@formula(event ~ 0 + fin+age+race+wexp+mar+paro+prio) ,rossi,
tol = 1e-8)

using BenchmarkTools # Expect around 4 ms (Matlab with Efron method is around 20)
@benchmark Survival.coxph(@formula(event ~ 0 + fin+age+race+wexp+mar+paro+prio) ,rossi,
tol = 1e-8)

@code_warntype coxph(@formula(event ~ 0 + fin+age+race+wexp+mar+paro+prio) ,rossi,
tol = 1e-8)



using Base.Test
using JLD
creload("Survival")

filepath = joinpath(Pkg.dir("Survival", "test"), "rossi.csv")
rossi = readtable(filepath)
rossi[:event] = Survival.EventTime.(rossi[:week],rossi[:arrest] .== 1)

outcome = Survival.coxph( @formula(event ~ 0+ fin+age+race+wexp+mar+paro+prio) ,rossi ; tol = 1e-8)
outcome_coefmat = coeftable(outcome)

filepath_coefs = joinpath(Pkg.dir("Survival", "test"), "expected_coefmat.jld")
expected_coefmat = JLD.load(filepath_coefs, "expected_coefmat")
@test outcome_coefmat.cols[1:3] ≈ expected_coefmat.cols[1:3] atol = 1e-6
