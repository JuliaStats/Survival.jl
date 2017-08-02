using Survival
using ClobberingReload
using DataFrames
creload("Survival")
filepath = joinpath(Pkg.dir("Survival", "test"), "rossi.csv")
rossi = readtable(filepath)
rossi[:event] = EventTime.(rossi[:week],rossi[:arrest] .== 1)

using BenchmarkTools # Expect around 4 ms (Matlab with Efron method is around 20)
@benchmark fit( Survival.CoxModel,@formula(event ~ fin+age+race+wexp+mar+paro+prio) ,rossi ; tol = 1e-8)


s = StatsBase.fit(Survival.CoxModel, @formula(event ~ 0 + fin+age+race+wexp+mar+paro+prio) ,rossi,
tol = 1e-8)

coxph(@formula(event ~ 0 + fin+age+race+wexp+mar+paro+prio) ,rossi,
tol = 1e-8)

coeftable(s).cols[1:3]

DataFrames.coefnames(s.mf)
s
s.mf.terms.intercept = false
coef(s)

mf = ModelFrame(f, df, contrasts=contrasts)
mm = ModelMatrix(mf)
y = model_response(mf)
