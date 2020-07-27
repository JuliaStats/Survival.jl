"""
    fit(CumulativeIncidence, CompetingEventTimes)
    fit(CumulativeIncidence, times, status, eventofinterest, censoringevent)

Given a vector of times to events, a vector of indicators for the event of interest,
and a vector of indicators for the competing event, compute the an estimate of the
cumulative incidence of failure.

The resulting `CumulativeIncidence` object has the following fields:

* `times`: Distinct event times
* `neventsofinterest`: Number of observed events of interest at each time
* `nallevents`: Number of observed events of any kind at each time
* `ncensor`: Number of right censorings for all events at each time
* `natrisk`: Size of the risk set at each time
* `estimate`: Estimate of the cumulative incidence of failure
* `stderr`: Standard error of the estimate of the cumulative incidence of failure
* `survival`: Estimate of the overall Kaplan-Meier survival
* `survivalstderr`: Standard error of the estimate of the overall Kaplan-Meier survival

### Formulas

In the presence of competing risks, `1-KM`, where `KM` is the Kaplan-Meier estimate,
is uninterpretable and is a biased estimate of the failure probability. The
cumulative incidence estimator of Kalbfleisch and Prentice (1980) is a function
of the hazards of both the event of interest and the competing event, and provides
an unbiased estimate of the failure probability.

The estimator of the cumulative incidence for event ``k`` is given by
``
\\hat{I}_k(t) = \\sum_{i: t_i < t} \\hat{S}(t_{i-1}) \\frac{d_{k,i}}{n_i}
``
where ``d_{k,i}`` are the events of interest ``k`` at time ``t_i``,
``n_i`` are the individuals at risk at that time,
and ``\\hat{S}(t_{i-1})`` is the usual Kaplan-Meier estimate of survival.

Standard errors are computed using the Delta method.

### References

* Kalbfleisch, J. D., and Prentice, R. L. (1980). *The Statistical Analysis of
  Failure Time Data*. New York, NY: John Wiley.
"""
struct CumulativeIncidence{T<:Real} <: NonparametricEstimator
    times::Vector{T}
    neventsofinterest::Vector{Int}
    nallevents::Vector{Int}
    ncensor::Vector{Int}
    natrisk::Vector{Int}
    estimate::Vector{Float64}
    stderr::Vector{Float64}
    survival::Vector{Float64} # overall KM estimator (needed for CIF, so might as well return it)
    survivalstderr::Vector{Float64} # se of overall KM estimator
end




estimator_start(::Type{CumulativeIncidence}) = 0.0  # Estimator starting point

estimator_update(::Type{CumulativeIncidence}, es, d, n, S) = es + S * (d / n) # Estimator update rule

stderr_start(::Type{CumulativeIncidence}) = 0.0


# Delta method variance formula (as opposed to Aalen counting process method)
function stderr_update(::Type{CumulativeIncidence}, es, dᵢ, dₐ, nₐ, surv, gw, vartmp1, vartmp2)
    nextvartmp1 = vartmp1 + es * dₐ / (nₐ * (nₐ - dₐ)) + surv * dᵢ / nₐ^2
    nextvartmp2 = vartmp2 + es^2 * dₐ / (nₐ * (nₐ - dₐ)) + surv^2 * (nₐ - dᵢ) * dᵢ / nₐ^3 + 2 * es * surv * dᵢ / nₐ^2
    vares = es^2 * gw - 2 * es * nextvartmp1 + nextvartmp2 
    return (vares,nextvartmp1,nextvartmp2)
end



# formula ref: Coviello and Boggess (2004), stcompet, The Stata Journal
# expand the variance formula into current values and sums of past values

# subscript "a" means any or all, while "i" means event of interest

function _estimator(::Type{CumulativeIncidence}, tte::AbstractVector{CompetingEventTime{T,S}}) where {T,S}
    nobs = length(tte)
    dᵢ = 0                   # Number of observed events of interest at time t
    dₐ = 0                   # Number of all observed events at time t
    cₐ = 0                   # Number of censored events at time t
    nₐ = nobs                # Number remaining at risk at time t
    es = estimator_start(CumulativeIncidence)  # Estimator starting point
    vares = stderr_start(CumulativeIncidence)     # Standard Error starting point
    surv = estimator_start(KaplanMeier)  # Survival Estimator starting point 
    gw = stderr_start(KaplanMeier)     # Standard Error starting point for Greenwood's portion

    times = T[]              # The set of unique event times
    neventsofinterest = Int[]  # Total observed events of interest at each time
    nallevents = Int[]          # Total observed events at each time
    ncensor = Int[]          # Total censored events at each time
    natrisk = Int[]          # Number at risk at each time
    estimator = Float64[]    # Estimates
    stderr = Float64[]       # Pointwise standard errors
    survival = Float64[]       # Overall survival function estimate
    survivalstderr = Float64[]       # Pointwise standard error of Overall survival function estimate

    t_prev = zero(T)
    vartmp1 = 0.0
    vartmp2 = 0.0

    @inbounds for i = 1:nobs
        tte_i = tte[i]
        t = tte_i.time
        s = iseventofinterest(tte_i)
        a = isevent(tte_i)
        # Aggregate over tied times
        if t == t_prev
            dᵢ += s
            cₐ += !a
            dₐ += a
            continue
        elseif !iszero(t_prev)
            gw = stderr_update(KaplanMeier, gw, dₐ, nₐ)
            es = estimator_update(CumulativeIncidence, es, dᵢ, nₐ, surv)
            vares,vartmp1,vartmp2 = stderr_update(CumulativeIncidence, es, dᵢ, dₐ, nₐ, surv, gw, vartmp1, vartmp2)
            surv = estimator_update(KaplanMeier, surv, dₐ, nₐ)
            
            push!(times, t_prev)
            push!(neventsofinterest, dᵢ)
            push!(nallevents, dₐ)
            push!(ncensor, cₐ)
            push!(natrisk, nₐ)
            push!(estimator, es)
            push!(stderr, sqrt(vares))
            push!(survival, surv)
            push!(survivalstderr, surv*sqrt(gw))
        end
        nₐ -= dₐ + cₐ
        cₐ = !a
        dₐ = a
        dᵢ = s
        t_prev = t
    end

    # We need to do this one more time to capture the last time
    # since everything in the loop is lagged
    push!(times, t_prev)
    push!(neventsofinterest, dᵢ)
    push!(nallevents, dₐ)
    push!(ncensor, cₐ)
    push!(natrisk, nₐ)
    push!(estimator, es)
    push!(stderr, sqrt(vares))
    push!(survival, surv)
    push!(survivalstderr, surv*sqrt(gw))

    return CumulativeIncidence{T}(times, neventsofinterest, nallevents, ncensor, natrisk, estimator, stderr, survival, survivalstderr)
end


function StatsBase.fit(::Type{CumulativeIncidence},tte::AbstractVector{CompetingEventTime{T,S}}) where {T<:Real,S}
    nobs = length(tte)
    if nobs == 0
        throw(ArgumentError("The sample must be nonempty."))
    end

    sortedevents = sort(tte)

    return _estimator(CumulativeIncidence, sortedevents)
end




# helper function for people too lazy to create CompetingEventTime objects
function StatsBase.fit(::Type{CumulativeIncidence},
                       times::AbstractVector{T},
                       status::AbstractVector{S},
                       eventofinterest::S,censoringevent::S) where {T<:Real,S}
    nobs = length(times)
    if !(nobs == length(status))
        throw(DimensionMismatch("The input vectors must have the same length"))
    end
    if nobs == 0
        throw(ArgumentError("the sample must be nonempty"))
    end
    tte = CompetingEventTime.(times,status,eventofinterest,censoringevent)
    return StatsBase.fit(CumulativeIncidence,tte)
end


"""
    confint(km::CumulativeIncidence, α=0.05)

Compute the pointwise log-log transformed confidence intervals for the 
cumulative incidence function as a vector of tuples.
"""
function StatsBase.confint(km::CumulativeIncidence, α::Float64=0.05)
    q = quantile(Normal(), 1 - α/2)
    return map(km.estimate, km.stderr) do es, se
        a = q * se / (es *log(es))
        es^exp(-a), es^exp(a)
    end
end




