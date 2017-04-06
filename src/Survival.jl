__precompile__()

module Survival

using StatsBase

export
    EventTime,
    EventTimeVector,
    isevent,
    iscensored,

    KaplanMeier,
    fit

include("eventtimes.jl")
include("kaplanmeier.jl")

end # module
