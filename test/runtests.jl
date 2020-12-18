using Survival
using Test
using CSV
using DataFrames
using Distributions
using LinearAlgebra
using StatsBase
using StatsModels

@testset "Event times" begin
    @test isevent(EventTime{Int}(44, true))
    @test !isevent(EventTime(3.2, false))
    @test isevent(EventTime(2.5f0))

    @test !iscensored(EventTime(3))
    @test iscensored(EventTime(2.1, false))

    @test eltype(EventTime(1)) == Int

    @test sprint(show, EventTime(1, true)) == "1"
    @test sprint(show, EventTime(1, false)) == "1+"

    @test convert(Int, EventTime(1)) == 1
    @test convert(Int, EventTime(1, false)) == 1
    @test convert(Float64, EventTime(1)) == 1.0
    @test convert(EventTime, 1) == EventTime(1)

    @test isless(EventTime(1), EventTime(1, false))
    @test !isless(EventTime(2), EventTime(1, false))

    let x = [EventTime(2, false), EventTime(1), EventTime(2)]
        @test sort(x) == [EventTime(1), EventTime(2), EventTime(2, false)]
    end
end

@testset "Kaplan-Meier" begin
    t = [
        310, 361, 654, 728,  61,  81, 520, 473, 107, 122, 965, 731, 153, 433, 145,  95, 765,
        735,   5, 687, 345, 444,  60, 208, 821, 305, 226, 426, 705, 363, 167, 641, 740, 245,
        588, 166, 559, 450, 529, 351, 201, 524, 199, 550, 551, 543, 293, 511, 511, 371, 201,
         62, 356, 340, 315, 182, 364, 376, 384, 268, 266, 194, 348, 382, 296, 186, 145, 269,
        350, 272, 292, 332, 285, 243, 276,  79, 240, 202, 235, 224, 239, 173, 252,  92, 192,
        211, 175, 203, 105, 177,
    ]
    s = [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1,
        1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1,
        0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
    ]

    # Results generated from R 3.3.2 using survival 2.41
    # survfit(Surv(t, s) ~ 1, conf.type = "log-log")
    # $surv
    r_surv = [
        0.9888889, 0.9777778, 0.9666667, 0.9555556, 0.9444444,  0.9333333,  0.9333333,  0.9220884,
        0.9220884, 0.9107045, 0.8993207, 0.8765531, 0.8651693,  0.8537855,  0.8424017,  0.8424017,
        0.8424017, 0.8424017, 0.8305369, 0.8186721, 0.8186721,  0.8066328,  0.7945935,  0.7705149,
        0.7705149, 0.7705149, 0.7580872, 0.7580872, 0.7580872,  0.7452383,  0.7452383,  0.7321639,
        0.7321639, 0.7321639, 0.7186054, 0.7186054, 0.7186054,  0.7045151,  0.7045151,  0.7045151,
        0.7045151, 0.6895254, 0.6895254, 0.6742026, 0.6742026,  0.6585235,  0.6428443,  0.6428443,
        0.6428443, 0.6263611, 0.609878,  0.5933948, 0.5769116,  0.5604284,  0.5604284,  0.5434457,
        0.526463,  0.526463,  0.5089143, 0.5089143, 0.5089143,  0.5089143,  0.4893406,  0.469767,
        0.4501934, 0.4306198, 0.4110461, 0.4110461, 0.3894121,  0.3677781,  0.3677781,  0.3677781,
        0.3432596, 0.3432596, 0.3432596, 0.3432596, 0.3120542,  0.2808487,  0.2496433,  0.2184379,
        0.1872325, 0.1560271, 0.1248217, 0.1248217, 0.08321444, 0.08321444, 0.08321444,
    ]
    # $n.risk
    r_risk = [
        90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68,
        67, 66, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45,
        44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23,
        22, 21, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,  9,  8,  7,  6,  5,  4,  3,  2,  1
    ]
    # $std.err
    r_stderr = [
        0.01117336, 0.01589104, 0.01957401, 0.02273314, 0.02556550, 0.02817181, 0.02817181,
        0.03066888, 0.03066888, 0.03308929, 0.03539956, 0.03977328, 0.04186640, 0.04391166,
        0.04591747, 0.04591747, 0.04591747, 0.04591747, 0.04805852, 0.05016633, 0.05016633,
        0.05230824, 0.05442696, 0.05861552, 0.05861552, 0.05861552, 0.06082918, 0.06082918,
        0.06082918, 0.06318557, 0.06318557, 0.06561783, 0.06561783, 0.06561783, 0.06822833,
        0.06822833, 0.06822833, 0.07104408, 0.07104408, 0.07104408, 0.07104408, 0.07422800,
        0.07422800, 0.07755545, 0.07755545, 0.08104663, 0.08455340, 0.08455340, 0.08455340,
        0.08845361, 0.09238657, 0.09636405, 0.10039761, 0.10449888, 0.10449888, 0.10893570,
        0.11346828, 0.11346828, 0.11842498, 0.11842498, 0.11842498, 0.11842498, 0.12475150,
        0.13126159, 0.13798985, 0.14497408, 0.15225631, 0.15225631, 0.16157339, 0.17138826,
        0.17138826, 0.17138826, 0.18475887, 0.18475887, 0.18475887, 0.18475887, 0.20791044,
        0.23310483, 0.26120251, 0.29340057, 0.33150176, 0.37845310, 0.43957565, 0.43957565,
        0.59991117, 0.59991117, 0.59991117,
    ]
    # $lower
    r_lower = [
        0.92374348, 0.91406011, 0.90021666, 0.88590944, 0.87167275, 0.85762000, 0.85762000,
        0.84350521, 0.84350521, 0.82935113, 0.81542216, 0.78813835, 0.77474465, 0.76149619,
        0.74838092, 0.74838092, 0.74838092, 0.74838092, 0.73463920, 0.72104592, 0.72104592,
        0.70732360, 0.69373775, 0.66693587, 0.66693587, 0.66693587, 0.65314492, 0.65314492,
        0.65314492, 0.63887126, 0.63887126, 0.62441211, 0.62441211, 0.62441211, 0.60940174,
        0.60940174, 0.60940174, 0.59378857, 0.59378857, 0.59378857, 0.59378857, 0.57705901,
        0.57705901, 0.56006894, 0.56006894, 0.54279401, 0.52574922, 0.52574922, 0.52574922,
        0.50779307, 0.49008837, 0.47261930, 0.45537265, 0.43833741, 0.43833741, 0.42084761,
        0.40357982, 0.40357982, 0.38579388, 0.38579388, 0.38579388, 0.38579388, 0.36559209,
        0.34575866, 0.32627779, 0.30713724, 0.28832788, 0.28832788, 0.26728290, 0.24672627,
        0.24672627, 0.24672627, 0.22307271, 0.22307271, 0.22307271, 0.22307271, 0.19157693,
        0.16205660, 0.13440887, 0.10859921, 0.08465965, 0.06269905, 0.04293072, 0.04293072,
        0.01850525, 0.01850525, 0.01850525,
    ]
    # $upper
    r_upper = [
        0.9984273, 0.9943955, 0.9891262, 0.9830833, 0.9764926, 0.9694844, 0.9694844,
        0.9620777, 0.9620777, 0.9543175, 0.9463099, 0.9296781, 0.9211003, 0.9123712,
        0.9035043, 0.9035043, 0.9035043, 0.9035043, 0.8942180, 0.8848013, 0.8848013,
        0.8751547, 0.8653903, 0.8455370, 0.8455370, 0.8455370, 0.8352056, 0.8352056,
        0.8352056, 0.8244966, 0.8244966, 0.8135325, 0.8135325, 0.8135325, 0.8021438,
        0.8021438, 0.8021438, 0.7902942, 0.7902942, 0.7902942, 0.7902942, 0.7777438,
        0.7777438, 0.7648356, 0.7648356, 0.7515503, 0.7381102, 0.7381102, 0.7381102,
        0.7240036, 0.7097282, 0.6952904, 0.6806956, 0.6659482, 0.6659482, 0.6507130,
        0.6353162, 0.6353162, 0.6193746, 0.6193746, 0.6193746, 0.6193746, 0.6019235,
        0.5842205, 0.5662707, 0.5480775, 0.5296427, 0.5296427, 0.5096032, 0.4892180,
        0.4892180, 0.4892180, 0.4666889, 0.4666889, 0.4666889, 0.4666889, 0.4401063,
        0.4122125, 0.3830497, 0.3526092, 0.3208327, 0.2876052, 0.2527363, 0.2527363,
        0.2123639, 0.2123639, 0.2123639,
    ]

    km = fit(KaplanMeier, t, s)
    jl_surv = km.survival

    @test length(jl_surv) == length(r_surv)
    @test jl_surv ≈ r_surv atol=1e-6

    @test km.times == sort!(unique(t))
    @test km.natrisk == r_risk
    @test km.nevents == [sum(s[t .== tᵢ]) for tᵢ in sort!(unique(t))]
    @test km.ncensor == [sum(iszero, s[t .== tᵢ]) for tᵢ in sort!(unique(t))]
    @test km.stderr ≈ r_stderr atol=1e-6

    conf = confint(km)
    jl_lower = first.(conf)
    jl_upper = last.(conf)

    @test jl_lower ≈ r_lower atol=1e-6
    @test jl_upper ≈ r_upper atol=1e-6

    @test_throws DimensionMismatch fit(KaplanMeier, [1, 2], [true])
    @test_throws ArgumentError fit(KaplanMeier, Float64[], Bool[])

    km_et = fit(KaplanMeier, EventTime.(t, Bool.(s)))
    @test all(f->getfield(km, f) ≈ getfield(km_et, f), fieldnames(KaplanMeier))

    @test_throws ArgumentError fit(KaplanMeier, EventTime{Int}[])
end

@testset "Nelson-Aalen" begin
    t = [
        310, 361, 654, 728,  61,  81, 520, 473, 107, 122, 965, 731, 153, 433, 145,  95, 765,
        735,   5, 687, 345, 444,  60, 208, 821, 305, 226, 426, 705, 363, 167, 641, 740, 245,
        588, 166, 559, 450, 529, 351, 201, 524, 199, 550, 551, 543, 293, 511, 511, 371, 201,
         62, 356, 340, 315, 182, 364, 376, 384, 268, 266, 194, 348, 382, 296, 186, 145, 269,
        350, 272, 292, 332, 285, 243, 276,  79, 240, 202, 235, 224, 239, 173, 252,  92, 192,
        211, 175, 203, 105, 177,
    ]
    s = [
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1,
        1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1,
        0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0,
    ]

    na = fit(NelsonAalen, t, s)
    km = fit(KaplanMeier, t, s)

    @test na.times == km.times
    @test na.nevents == km.nevents
    @test na.ncensor == km.ncensor
    @test na.natrisk == km.natrisk
    @test exp.(-na.chaz[1:50]) ≈ km.survival[1:50] rtol=1e-2
    @test na.stderr[1:50] ≈ km.stderr[1:50] rtol=2e-2
    na_conf = confint(na)
    na_lower, na_upper = getindex.(na_conf, 1), getindex.(na_conf, 2)
    @test cdf.(Normal.(na.chaz, na.stderr), na_lower) ≈ fill(0.025, length(na.chaz)) rtol=1e-8
    @test cdf.(Normal.(na.chaz, na.stderr), na_upper) ≈ fill(0.975, length(na.chaz)) rtol=1e-8
    na_conf = confint(na, 0.01)
    na_lower, na_upper = getindex.(na_conf, 1), getindex.(na_conf, 2)
    @test cdf.(Normal.(na.chaz, na.stderr), na_lower) ≈ fill(0.005, length(na.chaz)) rtol=1e-8
    @test cdf.(Normal.(na.chaz, na.stderr), na_upper) ≈ fill(0.995, length(na.chaz)) rtol=1e-8
end


@testset "Cox" begin
    rossi = CSV.read(joinpath(@__DIR__, "data", "rossi.csv"), DataFrame)
    rossi.event = EventTime.(rossi.week, rossi.arrest .== 1)

    outcome = coxph(@formula(event ~ fin + age + race + wexp + mar + paro + prio), rossi; tol=1e-8)
    outcome_coefmat = coeftable(outcome)

    regressor_matrix = Matrix(rossi[!, [:fin, :age, :race, :wexp, :mar, :paro, :prio]])
    event_vector = rossi.event

    outcome_without_formula = coxph(regressor_matrix, event_vector)

    @test sprint(show, outcome_without_formula) == chomp("""
        CoxModel{Float64}

        Coefficients:
        ──────────────────────────────────────────────
              Estimate  Std.Error    z value  Pr(>|z|)
        ──────────────────────────────────────────────
        x1  -0.379416   0.191379   -1.98253     0.0474
        x2  -0.0574299  0.0219988  -2.61059     0.0090
        x3   0.31392    0.307995    1.01924     0.3081
        x4  -0.14981    0.212226   -0.705898    0.4803
        x5  -0.433724   0.38187    -1.13579     0.2560
        x6  -0.0848615  0.195756   -0.433505    0.6646
        x7   0.091521   0.0286469   3.1948      0.0014
        ──────────────────────────────────────────────
        """)

    coef_matrix = ModelMatrix(ModelFrame(@formula(event ~ 0 + fin + age + race + wexp + mar + paro + prio), rossi)).m
    outcome_from_matrix     = coxph(coef_matrix, rossi.event; tol=1e-8, l2_cost=0)
    outcome_from_matrix32   = coxph(Float32.(coef_matrix), rossi.event; tol=1e-5)
    outcome_from_matrix_int = coxph(Int64.(coef_matrix), rossi.event; tol=1e-6, l2_cost=0.0)

    expected_coefs = [
        -0.379422   0.191379   -1.98256   0.0474;
        -0.0574377  0.0219995  -2.61087   0.0090;
         0.3139     0.307993    1.01918   0.3081;
        -0.149796   0.212224   -0.705837  0.4803;
        -0.433704   0.381868   -1.13574   0.2561;
        -0.0848711  0.195757   -0.433554  0.6646;
         0.0914971  0.0286485   3.19378   0.0014
    ]

    @test coef(outcome_from_matrix) ≈ coef(outcome) atol=1e-5
    @test coef(outcome_from_matrix) ≈ coef(outcome_from_matrix32) atol=1e-4
    @test coef(outcome_from_matrix) ≈ coef(outcome_from_matrix_int) atol=1e-5
    @test nobs(outcome) == size(rossi, 1)
    @test dof(outcome) == 7
    @test loglikelihood(outcome) > nullloglikelihood(outcome)
    @test all(x->x > 0, eigen(outcome.model.fischer_info).values)
    @test outcome.model.fischer_info * vcov(outcome) ≈ I atol=1e-10
    @test norm(outcome.model.score) < 1e-5
    @test hcat(outcome_coefmat.cols[1:3]...) ≈ expected_coefs[:,1:3] atol=1e-5

    outcome_fin = coxph(@formula(event ~ fin), rossi; tol=1e-8)
    @test coeftable(outcome_fin).rownms == ["fin"]
    outcome_finrace = coxph(@formula(event ~ fin * race), rossi; tol=1e-8)
    @test coeftable(outcome_finrace).rownms == ["fin", "race","fin & race"]
    transform!(rossi, :fin => categorical, renamecols = false)
    outcome_fincat = coxph(@formula(event ~ fin), rossi; tol=1e-8)
    @test coeftable(outcome_fincat).rownms == ["fin: 1"]
    @test coef(outcome_fin) ≈ coef(outcome_fincat) atol=1e-8
    outcome_fincatrace = coxph(@formula(event ~ fin * race), rossi; tol=1e-8)
    @test coeftable(outcome_fincatrace).rownms == ["fin: 1", "race","fin: 1 & race"]
    @test coef(outcome_fincatrace) ≈ coef(outcome_finrace) atol=1e-8
    transform!(rossi, :race => categorical, renamecols = false)
    outcome_fincatracecat = coxph(@formula(event ~ fin * race), rossi; tol=1e-8)
    @test coeftable(outcome_fincatracecat).rownms == ["fin: 1", "race: 1","fin: 1 & race: 1"]
    @test coef(outcome_fincatracecat) ≈ coef(outcome_finrace) atol=1e-8
end

@testset "Newton-Raphson" begin
    function fgh!(x, grad, hes, compute_ders)
        if compute_ders
            grad[1] = 2exp(x[1])*(exp(x[1]) - 1)
            hes[1,1] = 2exp(x[1])*(2exp(x[1]) - 1)
        end
        (exp(x[1]) - 1)^2
    end
    x, y, grad, hes = Survival.newton_raphson(fgh!, [2.2], tol=1e-5)
    @test x ≈ [0] atol=1e-5
    @test y ≈ 0 atol=1e-5
    @test grad ≈ [0] atol=1e-5
    @test hes ≈ [2] atol=1e-5
    @test_throws ConvergenceException Survival.newton_raphson(fgh!, [2.2], max_iter=2)

    function wrong_fgh!(x, grad, hes, compute_ders)
        if compute_ders
            grad[1] = -2exp(x[1])*(exp(x[1]) - 1) # wrong sign
            hes[1,1] = 2exp(x[1])*(2exp(x[1]) - 1)
        end
        (exp(x[1]) - 1)^2
    end
    @test_throws ErrorException Survival.newton_raphson(wrong_fgh!, [2.2])
end
