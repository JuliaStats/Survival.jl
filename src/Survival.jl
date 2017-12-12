__precompile__()

module Survival

using StatsBase, StatsModels
using Distributions
using PositiveFactorizations

export
    EventTime,
    isevent,
    iscensored,

    KaplanMeier,
    fit,
    confint,

    coxph,
    coef,
    loglikelihood,
    nullloglikelihood,
    nobs,
    dof,
    vcov,
    stderr


abstract type AbstractEstimator end
abstract type NonparametricEstimator <: AbstractEstimator end

include("eventtimes.jl")
include("kaplanmeier.jl")
include("cox.jl")
include("optimization.jl")

end # module
