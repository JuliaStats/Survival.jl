module Survival

using Compat
using Distributions
using LinearAlgebra
using Optim
using StatsAPI
using StatsBase
using StatsModels

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
    modelmatrix,
    loglikelihood,
    nullloglikelihood,
    nobs,
    dof,
    dof_residual,
    vcov,
    stderror


abstract type AbstractEstimator end
abstract type NonparametricEstimator <: AbstractEstimator end

include("eventtimes.jl")
include("estimator.jl")
include("kaplanmeier.jl")
include("nelsonaalen.jl")
include("cox.jl")

end # module
