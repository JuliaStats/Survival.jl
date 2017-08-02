filepath = joinpath(Pkg.dir("Survival", "test"), "rossi.csv")
rossi = readtable(filepath)
rossi[:event] = EventTime.(rossi[:week],rossi[:arrest] .== 1)

outcome = coxph( @formula(event ~ 0+ fin+age+race+wexp+mar+paro+prio) ,rossi ; tol = 1e-8)
outcome_coefmat = coeftable(outcome)

filepath_coefs = joinpath(Pkg.dir("Survival", "test"), "expected_coefmat.jld")
expected_coefmat = JLD.load(filepath_coefs, "expected_coefmat")
@test outcome_coefmat.cols[1:3] ≈ expected_coefmat.cols[1:3] atol = 1e-6
println("Cox regression is fine")
println("$(outcome_coefmat) \n ≈ \n $expected_coefmat")
