using Survival
using ClobberingReload
using DataFrames
creload("Survival")
filepath = joinpath(Pkg.dir("Survival", "test"), "rossi.csv")
rossi = readtable(filepath)
rossi[:event] = EventTime.(rossi[:week],rossi[:arrest] .== 1)
Survival.coxph( @formula(event ~ fin+age+race+wexp+mar+paro+prio) ,rossi ; tol = 1e-8)
