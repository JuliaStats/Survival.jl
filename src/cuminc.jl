"""
    fit(CumulativeIncidence, times, events, competing)

Given a vector of times to events, a vector of indicators for the event of interest,
and a vector of indicators for the competing event, compute the an estimate of the
cumulative incidence of failure.

The resulting `CumulativeIncidence` object has the following fields:

* `times`: Distinct event times
* `nevents`: Number of observed events of interest at each time
* `ncensor`: Number of right censorings for the event of interest at each time
* `ncompete`: Number of observed competing events at each time
* `natrisk`: Size of the risk set at each time
* `cuminc`: Estimate of the cumulative incidence of failure

### Formulas

In the presence of competing risks, `1-KM`, where `KM` is the Kaplan-Meier estimate,
is uninterpretable and is a biased estimate of the failure probability. The
cumulative incidence estimator of Kalbfleisch and Prentice (1980) is a function
of the hazards of both the event of interest and the competing event, and provides
an unbiased estimate of the failure probability.

The estimator is given by
``
\\hat{I}(t) = \\sum_{i: t_i < t} \\frac{d_i}{n_i} \\hat{S}_1(t) \\hat{S}_2(t)
``
where ``\\hat{S}_1(t)`` is the Kaplan-Meier estimator using the event of interest
and ``\\hat{S}_2(t)`` is the K-M estimator using the competing event.

### References

* Kalbfleisch, J. D., and Prentice, R. L. (1980). *The Statistical Analysis of
  Failure Time Data*. New York, NY: John Wiley.
"""
struct CumulativeIncidence{T<:Real} <: NonparametricEstimator
    times::Vector{T}
    nevents::Vector{Int}
    ncensor::Vector{Int}
    ncompete::Vector{Int}
    natrisk::Vector{Int}
    cuminc::Vector{Float64}
end

function _cuminc(tte::AbstractVector{T}, status::BitVector, compete::BitVector) where {T}
    nobs = length(tte)
    dᵢ = 0              # Number of observed events of interest at time t
    cᵢ = 0              # Number of censored events of interest at time t
    cr = 0              # Number censored for the competing event
    rᵢ = 0              # Number of observed competing events at time t
    nᵢ = nobs           # Number remaining at risk just before time t
    nr = nobs           # Number remaining at risk for the competing event
    km1 = 1.0           # Ŝ(t) for the event of interest
    km2 = 1.0           # Ŝ(t) for the competing event
    inc = 0.0           # Î(t)

    times = T[]         # The set of unique event times
    nevents = Int[]     # Total observed events of interest at each time
    ncensor = Int[]     # Total censored events of interest at each time
    ncompete = Int[]    # Total observed competing events at each time
    natrisk = Int[]     # Number at risk at each time
    cuminc = Float64[]  # Cumulative incidence estimates

    t_prev = zero(T)

    @inbounds for i = 1:nobs
        t = tte[i]
        s = status[i]
        r = compete[i]
        # Aggregate over tied times
        if t == t_prev
            dᵢ += s
            cᵢ += !s
            rᵢ += r
            cr += !r
            continue
        elseif !iszero(t_prev)
            km1 *= 1 - dᵢ / nᵢ
            km2 *= 1 - rᵢ / nr
            inc += km1 * km2 * dᵢ / nᵢ
            push!(times, t_prev)
            push!(nevents, dᵢ)
            push!(ncensor, cᵢ)
            push!(ncompete, rᵢ)
            push!(natrisk, nᵢ)
            push!(cuminc, inc)
        end
        nᵢ -= dᵢ + cᵢ
        nr -= rᵢ + cr
        dᵢ = s
        cᵢ = !s
        rᵢ = r
        cr = !r
        t_prev = t
    end

    push!(times, t_prev)
    push!(nevents, dᵢ)
    push!(ncensor, cᵢ)
    push!(ncompete, rᵢ)
    push!(natrisk, nᵢ)
    push!(cuminc, inc)

    return CumulativeIncidence{T}(times, nevents, ncensor, ncompete, natrisk, cuminc)
end

function StatsBase.fit(::Type{CumulativeIncidence},
                       times::AbstractVector{T},
                       status::AbstractVector{<:Integer},
                       compete::AbstractVector{<:Integer}) where {T<:Real}
    nobs = length(times)
    if !(nobs == length(status) == length(compete))
        throw(DimensionMismatch("the input vectors must have the same length"))
    end
    if nobs == 0
        throw(ArgumentError("the sample must be nonempty"))
    end
    p = sortperm(times)
    t = times[p]
    s = BitVector(status[p])
    r = BitVector(compete[p])
    return _cuminc(t, s, r)
end
