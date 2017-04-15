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
    "text": "This package provides types and methods for performing survival analysis in Julia.Pages = [\n    \"km.md\",\n]\nDepth = 1"
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
    "text": "The Kaplan-Meier estimator is a nonparametric estimator of the survivor function, i.e. the probability of survival beyond a given time.The estimate is given by hatS(t) = prod_i t_i  t left( 1 - fracd_in_i right) where d_i is the number of observed events at time t_i and n_i is the number of subjects at risk for the event just before time t_i.The pointwise standard error of the log of the survivor function can be computed using Greenwood's formula: textSE(log hatS(t)) = sqrtsum_i t_i  t fracd_in_i (n_i - d_i)"
},

{
    "location": "km.html#API-1",
    "page": "Kaplan-Meier",
    "title": "API",
    "category": "section",
    "text": "Survival.KaplanMeier\nSurvival.fit(::Type{Survival.KaplanMeier}, ::Any, ::Any)\nSurvival.confint(::Survival.KaplanMeier, ::Any)"
},

{
    "location": "km.html#References-1",
    "page": "Kaplan-Meier",
    "title": "References",
    "category": "section",
    "text": "Kaplan, E. L., and Meier, P. (1958). Nonparametric Estimation from Incomplete Observations. Journal of the American Statistical Association, 53(282), 457-481. doi:10.2307/2281868\nGreenwood, M. (1926). A Report on the Natural Duration of Cancer. Reports on Public Health and Medical Subjects. London: Her Majesty's Stationery Office. 33, 1-26."
},

]}
