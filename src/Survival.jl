__precompile__()

module Survival

using StatsBase, StatsModels
using Distributions
using PositiveFactorizations

export
    EventTime,
    isevent,
    iscensored,

    NonparametricEstimator,
    KaplanMeier,
    NelsonAalen,
    fit,
    confint,

    CoxModel,
    coxph,
    coef,
    loglikelihood,
    nullloglikelihood,
    nobs,
    dof,
    vcov,
    stderr,
    baseline_hazard,
    baseline_cumulative_hazard,
    baseline_survival


abstract type AbstractEstimator end
abstract type NonparametricEstimator <: AbstractEstimator end

include("eventtimes.jl")
include("estimator.jl")
include("kaplanmeier.jl")
include("nelsonaalen.jl")
include("cox.jl")
include("optimization.jl")

end # module
