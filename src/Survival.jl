module Survival

using Compat
using Distributions
using LinearAlgebra
using StatsAPI
using StatsBase
using StatsModels
using Tables

using NLSolversBase: NLSolversBase
using Optim: Optim

export
    EventTime,
    EventTable,
    isevent,
    iscensored,

    NonparametricEstimator,
    KaplanMeier,
    NelsonAalen,
    fit,
    confint,

    kaplanmeier,
    censorticks,
    risktable_at,
    add_risktable!,

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
    stderror,
    confint


abstract type AbstractEstimator end
abstract type NonparametricEstimator <: AbstractEstimator end

include("eventtimes.jl")
include("estimator.jl")
include("kaplanmeier.jl")
include("nelsonaalen.jl")
include("cox.jl")

## AlgebraOfGraphics extension stubs (methods added by SurvivalAlgebraOfGraphicsExt)
function kaplanmeier end
function censorticks end
function risktable_at end
function add_risktable! end

## Deprecations

@noinline function StatsAPI.confint(estimator::NonparametricEstimator, alpha::Float64)
    Base.depwarn("`confint(estimator, alpha)` is deprecated, use " *
                 "`confint(estimator; level=alpha)` instead.", :confint)
    return confint(estimator; level=alpha)
end

end # module
