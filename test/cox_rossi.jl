filepath = joinpath(Pkg.dir("Survival", "test"), "rossi.csv")
rossi = readtable(filepath)
rossi[:event] = EventTime.(rossi[:week],rossi[:arrest] .== 1)

outcome = Survival.coxph( @formula(event ~ fin+age+race+wexp+mar+paro+prio) ,rossi ; tol = 1e-8)

filepath_coefs = joinpath(Pkg.dir("Survival", "test"), "expected_coefmat.jld")
expected_coefmat = JLD.load(filepath_coefs, "expected_coefmat")
@test outcome.coefmat.cols[1:3] ≈ expected_coefmat.cols[1:3] atol = 1e-6
println("Cox regression is fine")
println("$(outcome.coefmat) \n ≈ \n $expected_coefmat")
