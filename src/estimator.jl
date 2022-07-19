# Estimating functions with the following assumptions:
#  * The input is nonempty
#  * Time 0 is not included

function _n_unique_times(ets)
    # We know that `ets` is nonempty and sorted ascending by the `time` fields of the
    # elements, so we can count the unique values in `O(length(ets))` time with no
    # allocations by only comparing to the last observed unique time
    t_prev = first(ets).time
    n = 1
    for et in Iterators.drop(ets, 1)
        t = et.time
        if t != t_prev
            n += 1
            t_prev = t
        end
    end
    return n
end

function _estimator(::Type{S}, ets::AbstractVector{EventTime{T}}) where {S,T}
    outlen = _n_unique_times(ets)

    nobs = length(ets)
    dᵢ::Int = 0              # Number of observed events at time t
    cᵢ::Int = 0              # Number of censored events at time t
    nᵢ::Int = nobs           # Number remaining at risk at time t
    es = estimator_start(S)  # Estimator starting point
    gw = stderr_start(S)     # Standard Error starting point

    times = Vector{T}(undef, outlen)            # The set of unique event times
    nevents = Vector{Int}(undef, outlen)        # Total observed events at each time
    ncensor = Vector{Int}(undef, outlen)        # Total censored events at each time
    natrisk = Vector{Int}(undef, outlen)        # Number at risk at each time
    estimator = Vector{Float64}(undef, outlen)  # Estimates
    stderr = Vector{Float64}(undef, outlen)     # Pointwise standard errors

    t_prev = zero(T)
    outind = 1

    @inbounds for obsind in eachindex(ets)
        t = ets[obsind].time
        s = ets[obsind].status
        # Aggregate over tied times
        if t == t_prev
            dᵢ += s
            cᵢ += !s
            continue
        elseif !iszero(t_prev)
            es = estimator_update(S, es, dᵢ, nᵢ)
            gw = stderr_update(S, gw, dᵢ, nᵢ)
            times[outind] = t_prev
            nevents[outind] = dᵢ
            ncensor[outind] = cᵢ
            natrisk[outind] = nᵢ
            estimator[outind] = es
            stderr[outind] = sqrt(gw)
            outind += 1
        end
        nᵢ -= dᵢ + cᵢ
        dᵢ = s
        cᵢ = !s
        t_prev = t
    end

    # We need to do this one more time to capture the last time
    # since everything in the loop is lagged
    times[outind] = t_prev
    nevents[outind] = dᵢ
    ncensor[outind] = cᵢ
    natrisk[outind] = nᵢ
    estimator[outind] = es
    stderr[outind] = sqrt(gw)

    return S{T}(times, nevents, ncensor, natrisk, estimator, stderr)
end

function StatsBase.fit(::Type{S},
                       times::AbstractVector{T},
                       status::AbstractVector{<:Integer}) where {S<:NonparametricEstimator,T}
    ntimes = length(times)
    nstatus = length(status)
    if ntimes != nstatus
        throw(DimensionMismatch("number of event statuses does not match number of " *
                                "event times; got $nstatus and $ntimes, respectively"))
    end
    return fit(S, map(EventTime{T}, times, status))
end

function StatsBase.fit(::Type{S}, ets::AbstractVector{<:EventTime}) where S<:NonparametricEstimator
    isempty(ets) && throw(ArgumentError("can't compute $(nameof(S)) from 0 observations"))
    return _estimator(S, issorted(ets) ? ets : sort(ets))
end
