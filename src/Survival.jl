__precompile__()

module Survival

using StatsBase
using Distributions

export
    EventTime,
    isevent,
    iscensored,

    KaplanMeier,
    fit,
    confint

abstract type AbstractEstimator end
abstract type NonparametricEstimator <: AbstractEstimator end

include("eventtimes.jl")
include("kaplanmeier.jl")

end # module
