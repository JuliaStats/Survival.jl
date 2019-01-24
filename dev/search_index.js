var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "DocTestSetup = :(using Survival, StatsBase)\nCurrentModule = Survival"
},

{
    "location": "#Survival.jl-1",
    "page": "Home",
    "title": "Survival.jl",
    "category": "section",
    "text": "This package provides types and methods for performing survival analysis in Julia."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is not yet registered in Julia\'s General package registry, and so it must be installed using Pkg.add(\"https://github.com/ararslan/Survival.jl\")."
},

{
    "location": "#Contents-1",
    "page": "Home",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\n    \"events.md\",\n    \"km.md\",\n    \"cox.md\",\n]\nDepth = 1"
},

{
    "location": "events/#",
    "page": "Event Times",
    "title": "Event Times",
    "category": "page",
    "text": ""
},

{
    "location": "events/#Survival.EventTime",
    "page": "Event Times",
    "title": "Survival.EventTime",
    "category": "type",
    "text": "EventTime{T}\n\nImmutable object containing the real-valued time to an event as well as an indicator of whether the time corresponds to an observed event (true) or right censoring (false).\n\n\n\n\n\n"
},

{
    "location": "events/#Event-Times-1",
    "page": "Event Times",
    "title": "Event Times",
    "category": "section",
    "text": "A crucial concept in survival analysis is the time elapsed between some landmark and a particular event of interest. As an example, say you\'re running a clinical trial to investigate the efficacy of a new anticonvulsant. You may be interested in the time from the start of therapy to the first epileptic seizure for each patient. But if a patient dies or otherwise goes off study before they have a seizure, you\'ll assume that a seizure would have occurred eventually, but you don\'t know when exactly. In this case the event time is right censored; the true event time is unknown, all you know is that it exceeds the observed time.A dedicated type is provided to conveniently store right censored data.Survival.EventTime"
},

{
    "location": "km/#",
    "page": "Kaplan-Meier",
    "title": "Kaplan-Meier",
    "category": "page",
    "text": ""
},

{
    "location": "km/#Kaplan-Meier-Estimator-1",
    "page": "Kaplan-Meier",
    "title": "Kaplan-Meier Estimator",
    "category": "section",
    "text": "The Kaplan-Meier estimator is a nonparametric estimator of the survivor function, i.e. the probability of survival beyond a given time.The estimate is given byhatS(t) = prod_i t_i  t left( 1 - fracd_in_i right)where d_i is the number of observed events at time t_i and n_i is the number of subjects at risk for the event just before time t_i.The pointwise standard error of the log of the survivor function can be computed using Greenwood\'s formula:textSE(log hatS(t)) = sqrtsum_i t_i  t fracd_in_i (n_i - d_i)"
},

{
    "location": "km/#Survival.KaplanMeier",
    "page": "Kaplan-Meier",
    "title": "Survival.KaplanMeier",
    "category": "type",
    "text": "KaplanMeier\n\nAn immutable type containing survivor function estimates computed using the Kaplan-Meier method. The type has the following fields:\n\ntimes: Distinct event times\nnevents: Number of observed events at each time\nncensor: Number of right censored events at each time\nnatrisk: Size of the risk set at each time\nsurvival: Estimate of the survival probability at each time\nstderr: Standard error of the log survivor function at each time\n\nUse fit(KaplanMeier, ...) to compute the estimates and construct this type.\n\n\n\n\n\n"
},

{
    "location": "km/#StatsBase.fit-Tuple{Type{KaplanMeier},Any,Any}",
    "page": "Kaplan-Meier",
    "title": "StatsBase.fit",
    "category": "method",
    "text": "fit(KaplanMeier, times, status) -> KaplanMeier\n\nGiven a vector of times to events and a corresponding vector of indicators that dictate whether each time is an observed event or is right censored, compute the Kaplan-Meier estimate of the survivor function.\n\n\n\n\n\n"
},

{
    "location": "km/#StatsBase.confint-Tuple{KaplanMeier,Float64}",
    "page": "Kaplan-Meier",
    "title": "StatsBase.confint",
    "category": "method",
    "text": "confint(km::KaplanMeier, α=0.05)\n\nCompute the pointwise log-log transformed confidence intervals for the survivor function as a vector of tuples.\n\n\n\n\n\n"
},

{
    "location": "km/#API-1",
    "page": "Kaplan-Meier",
    "title": "API",
    "category": "section",
    "text": "Survival.KaplanMeier\nStatsBase.fit(::Type{KaplanMeier}, ::Any, ::Any)\nStatsBase.confint(::KaplanMeier, ::Float64)"
},

{
    "location": "km/#References-1",
    "page": "Kaplan-Meier",
    "title": "References",
    "category": "section",
    "text": "Kaplan, E. L., and Meier, P. (1958). Nonparametric Estimation from Incomplete Observations. Journal of the American Statistical Association, 53(282), 457-481. doi:10.2307/2281868\nGreenwood, M. (1926). A Report on the Natural Duration of Cancer. Reports on Public Health and Medical Subjects. London: Her Majesty\'s Stationery Office. 33, 1-26."
},

{
    "location": "na/#",
    "page": "Nelson-Aalen",
    "title": "Nelson-Aalen",
    "category": "page",
    "text": ""
},

{
    "location": "na/#Nelson-Aalen-Estimator-1",
    "page": "Nelson-Aalen",
    "title": "Nelson-Aalen Estimator",
    "category": "section",
    "text": "The Nelson-Aalen estimator is a nonparametric estimator of the cumulative hazard function.The estimate is given byhatH(t) = sum_i t_i  t fracd_in_iwhere d_i is the number of observed events at time t_i and n_i is the number of subjects at risk for the event just before time t_i.The pointwise standard error of the log of the survivor function can be computed directly as the standard error or a Bernoulli random variable with d_i successes from n_i samples:textSE(hatH(t)) = sqrtsum_i t_i  t fracd_i(n_i-d_i)n_i^3"
},

{
    "location": "na/#Survival.NelsonAalen",
    "page": "Nelson-Aalen",
    "title": "Survival.NelsonAalen",
    "category": "type",
    "text": "NelsonAalen\n\nAn immutable type containing cumulative hazard function estimates computed using the Nelson-Aalen method. The type has the following fields:\n\ntimes: Distinct event times\nnevents: Number of observed events at each time\nncensor: Number of right censored events at each time\nnatrisk: Size of the risk set at each time\nchaz: Estimate of the cumulative hazard at each time\nstderr: Standard error of the cumulative hazard\n\nUse fit(NelsonAalen, ...) to compute the estimates and construct this type.\n\n\n\n\n\n"
},

{
    "location": "na/#StatsBase.fit-Tuple{Type{NelsonAalen},Any,Any}",
    "page": "Nelson-Aalen",
    "title": "StatsBase.fit",
    "category": "method",
    "text": "fit(NelsonAalen, times, status) -> NelsonAalen\n\nGiven a vector of times to events and a corresponding vector of indicators that dictate whether each time is an observed event or is right censored, compute the Nelson-Aalen estimate of the cumulative hazard rate function.\n\n\n\n\n\n"
},

{
    "location": "na/#StatsBase.confint-Tuple{NelsonAalen,Float64}",
    "page": "Nelson-Aalen",
    "title": "StatsBase.confint",
    "category": "method",
    "text": "confint(na::NelsonAalen, α=0.05)\n\nCompute the pointwise confidence intervals for the cumulative hazard function as a vector of tuples.\n\n\n\n\n\n"
},

{
    "location": "na/#API-1",
    "page": "Nelson-Aalen",
    "title": "API",
    "category": "section",
    "text": "Survival.NelsonAalen\nStatsBase.fit(::Type{NelsonAalen}, ::Any, ::Any)\nStatsBase.confint(::NelsonAalen, ::Float64)"
},

{
    "location": "na/#References-1",
    "page": "Nelson-Aalen",
    "title": "References",
    "category": "section",
    "text": "Nelson, W. (1969). Hazard plotting for incomplete failure data. Journal of Quality Technology 1, 27–52."
},

{
    "location": "cox/#",
    "page": "Cox",
    "title": "Cox",
    "category": "page",
    "text": ""
},

{
    "location": "cox/#Cox-Proportional-Hazards-Model-1",
    "page": "Cox",
    "title": "Cox Proportional Hazards Model",
    "category": "section",
    "text": "The Cox proportional hazards model is a semiparametric regression model used to fit survival models without knowing the distribution. It is based on the assumption that covariates affect the hazard function multiplicatively. That is,lambda(t  X_i) = lambda_0(t) exp(X_i cdot beta)where lambda(tX_i) is the estimated hazard for sample i, lambda_0 is the baseline hazard, X_i is the vector of covariates for sample i, and beta is the vector of coefficients in the model."
},

{
    "location": "cox/#StatsBase.fit-Tuple{Type{CoxModel},AbstractArray{T,2} where T,AbstractArray{T,1} where T}",
    "page": "Cox",
    "title": "StatsBase.fit",
    "category": "method",
    "text": "fit(::Type{CoxModel}, M::AbstractMatrix, y::AbstractVector; kwargs...)\n\nGiven a matrix M of predictors and a corresponding vector of events, compute the Cox proportional hazard model estimate of coefficients. Returns a CoxModel object.\n\n\n\n\n\n"
},

{
    "location": "cox/#API-1",
    "page": "Cox",
    "title": "API",
    "category": "section",
    "text": "StatsBase.fit(::Type{CoxModel}, M::AbstractMatrix, y::AbstractVector; kwargs...)"
},

{
    "location": "cox/#References-1",
    "page": "Cox",
    "title": "References",
    "category": "section",
    "text": "Cox, D. R. (1972). Regression models and life tables (with discussion). Journal of the Royal Statistical Society, Series B, 34:187–220."
},

]}
