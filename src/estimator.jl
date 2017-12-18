# Estimating functions with the following assumptions:
#  * The input is nonempty
#  * Time 0 is not included

function _estimator(::Type{S}, tte::AbstractVector{T}, status::BitVector) where {S, T}
    nobs = length(tte)
    dᵢ = 0                   # Number of observed events at time t
    cᵢ = 0                   # Number of censored events at time t
    nᵢ = nobs                # Number remaining at risk at time t
    es = estimator_start(S)  # Estimator starting point
    gw = stderr_start(S)     # Standard Error starting point

    times = T[]              # The set of unique event times
    nevents = Int[]          # Total observed events at each time
    ncensor = Int[]          # Total censored events at each time
    natrisk = Int[]          # Number at risk at each time
    estimator = Float64[]    # Estimates
    stderr = Float64[]       # Pointwise standard errors

    t_prev = zero(T)

    @inbounds for i = 1:nobs
        t = tte[i]
        s = status[i]
        # Aggregate over tied times
        if t == t_prev
            dᵢ += s
            cᵢ += !s
            continue
        elseif !iszero(t_prev)
            es = estimator_update(S, es, dᵢ, nᵢ)
            gw = stderr_update(S, gw, dᵢ, nᵢ)
            push!(times, t_prev)
            push!(nevents, dᵢ)
            push!(ncensor, cᵢ)
            push!(natrisk, nᵢ)
            push!(estimator, es)
            push!(stderr, sqrt(gw))
        end
        nᵢ -= dᵢ + cᵢ
        dᵢ = s
        cᵢ = !s
        t_prev = t
    end

    # We need to do this one more time to capture the last time
    # since everything in the loop is lagged
    push!(times, t_prev)
    push!(nevents, dᵢ)
    push!(ncensor, cᵢ)
    push!(natrisk, nᵢ)
    push!(estimator, es)
    push!(stderr, sqrt(gw))

    return S{T}(times, nevents, ncensor, natrisk, estimator, stderr)
end

"""
    fit(::Type{S}, times, status) where S<:NonparametricEstimator

Given a vector of times to events and a corresponding vector of indicators that
dictate whether each time is an observed event or is right censored, compute the
model of type `S`. Return an object of type `S`: [`KaplanMeier`](@ref) and 
[`NelsonAalen`](@ref) are supported so far.
"""
function StatsBase.fit(::Type{S},
                       times::AbstractVector{T},
                       status::AbstractVector{<:Integer}) where {S<:NonparametricEstimator, T}
    nobs = length(times)
    if length(status) != nobs
        throw(DimensionMismatch("there must be as many event statuses as times"))
    end
    if nobs == 0
        throw(ArgumentError("the sample must be nonempty"))
    end
    p = sortperm(times)
    t = times[p]
    s = BitVector(status[p])
    return _estimator(S, t, s)
end

function StatsBase.fit(::Type{S}, ets::AbstractVector{<:EventTime}) where S<:NonparametricEstimator
    length(ets) > 0 || throw(ArgumentError("the sample must be nonempty"))
    x = sort(ets)
    # TODO: Refactor, since iterating over the EventTime objects directly in
    # the _km loop may actually be easier/more efficient than working with
    # the times and statuses as separate vectors. Plus it might be nice to
    # make this method the One True Method™ so that folks are encouraged to
    # use EventTimes instead of raw values.
    return fit(S, map(t->t.time, x), BitVector(map(t->t.status, x)))
end