__precompile__()

module Survival

using StatsBase
using Distributions
using PositiveFactorizations

export
    EventTime,
    isevent,
    iscensored,

    KaplanMeier,
    fit,
    confint,

    coxph

abstract type AbstractEstimator end
abstract type NonparametricEstimator <: AbstractEstimator end

include("eventtimes.jl")
include("kaplanmeier.jl")
include("cox_utils.jl")
include("cox_model.jl")
include("coxph.jl")
include("optimization.jl")

end # module
