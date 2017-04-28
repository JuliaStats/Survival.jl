var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Survival.jl-1",
    "page": "Home",
    "title": "Survival.jl",
    "category": "section",
    "text": "This package provides types and methods for performing survival analysis in Julia."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is not yet registered in Julia's package registry, and so it must be installed usingPkg.clone(\"https://github.com/ararslan/Survival.jl.git\")"
},

{
    "location": "index.html#Contents-1",
    "page": "Home",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\n    \"km.md\",\n]\nDepth = 1"
},

{
    "location": "km.html#",
    "page": "Kaplan-Meier",
    "title": "Kaplan-Meier",
    "category": "page",
    "text": ""
},

{
    "location": "km.html#Kaplan-Meier-estimator-1",
    "page": "Kaplan-Meier",
    "title": "Kaplan-Meier estimator",
    "category": "section",
    "text": "The Kaplan-Meier estimator is a nonparametric estimator of the survivor function, i.e. the probability of survival beyond a given time.The estimate is given byhatS(t) = prod_i t_i  t left( 1 - fracd_in_i right)where d_i is the number of observed events at time t_i and n_i is the number of subjects at risk for the event just before time t_i.The pointwise standard error of the log of the survivor function can be computed using Greenwood's formula:textSE(log hatS(t)) = sqrtsum_i t_i  t fracd_in_i (n_i - d_i)"
},

{
    "location": "km.html#Survival.KaplanMeier",
    "page": "Kaplan-Meier",
    "title": "Survival.KaplanMeier",
    "category": "Type",
    "text": "KaplanMeier\n\nAn immutable type containing survivor function estimates computed using the Kaplan-Meier method. The type has the following fields:\n\ntimes: Distinct event times\nnevents: Number of observed events at each time\nncensor: Number of right censored events at each time\nnatrisk: Size of the risk set at each time\nsurvival: Estimate of the survival probability at each time\nstderr: Standard error of the log survivor function at each time\n\nUse fit(KaplanMeier, ...) to compute the estimates and construct this type.\n\n\n\n"
},

{
    "location": "km.html#StatsBase.fit",
    "page": "Kaplan-Meier",
    "title": "StatsBase.fit",
    "category": "Function",
    "text": "fit(::Type{KaplanMeier}, times, status)\n\nGiven a vector of times to events and a corresponding vector of indicators that dictate whether each time is an observed event or is right censored, compute the Kaplan-Meier estimate of the survivor function. Returns a KaplanMeier object.\n\n\n\n"
},

{
    "location": "km.html#StatsBase.confint",
    "page": "Kaplan-Meier",
    "title": "StatsBase.confint",
    "category": "Function",
    "text": "confint(km::KaplanMeier, α=0.05)\n\nCompute the pointwise log-log transformed confidence intervals for the survivor function as a vector of tuples.\n\n\n\n"
},

{
    "location": "km.html#API-1",
    "page": "Kaplan-Meier",
    "title": "API",
    "category": "section",
    "text": "Survival.KaplanMeier\nStatsBase.fit\nStatsBase.confint"
},

{
    "location": "km.html#References-1",
    "page": "Kaplan-Meier",
    "title": "References",
    "category": "section",
    "text": "Kaplan, E. L., and Meier, P. (1958). Nonparametric Estimation from Incomplete Observations. Journal of the American Statistical Association, 53(282), 457-481. doi:10.2307/2281868\nGreenwood, M. (1926). A Report on the Natural Duration of Cancer. Reports on Public Health and Medical Subjects. London: Her Majesty's Stationery Office. 33, 1-26."
},

]}
