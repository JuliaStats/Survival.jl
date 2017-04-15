"""
    KaplanMeier

An immutable type containing survivor function estimates computed
using the Kaplan-Meier method.
The type has the following fields:

* `times`: Distinct event times
* `nevents`: Number of observed events at each time
* `ncensor`: Number of right censored events at each time
* `natrisk`: Size of the risk set at each time
* `survival`: Estimate of the survival probability at each time
* `stderr`: Standard error of the log survivor function at each time

Use `fit(KaplanMeier, ...)` to compute the estimates and construct
this type.
"""
struct KaplanMeier{T<:Real} <: NonparametricEstimator
    times::Vector{T}
    nevents::Vector{Int}
    ncensor::Vector{Int}
    natrisk::Vector{Int}
    survival::Vector{Float64}
    stderr::Vector{Float64}
end

# Internal Kaplan-Meier function with the following assumptions:
#  * The input is nonempty
#  * Time 0 is not included
function _km(tte::AbstractVector{T}, status::BitVector) where {T}
    nobs = length(tte)
    dᵢ = 0                # Number of observed events at time t
    cᵢ = 0                # Number of censored events at time t
    nᵢ = nobs             # Number remaining at risk at time t
    km = 1.0              # Ŝ(t)
    gw = 0.0              # SE(log(Ŝ(t)))

    times = T[]           # The set of unique event times
    nevents = Int[]       # Total observed events at each time
    ncensor = Int[]       # Total censored events at each time
    natrisk = Int[]       # Number at risk at each time
    survival = Float64[]  # Survival estimates
    stderr = Float64[]    # Pointwise standard errors for log(Ŝ(t))

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
            km *= 1 - dᵢ / nᵢ
            gw += dᵢ / (nᵢ * (nᵢ - dᵢ))
            push!(times, t_prev)
            push!(nevents, dᵢ)
            push!(ncensor, cᵢ)
            push!(natrisk, nᵢ)
            push!(survival, km)
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
    push!(survival, km)
    push!(stderr, sqrt(gw))

    return KaplanMeier{T}(times, nevents, ncensor, natrisk, survival, stderr)
end

"""
    fit(::Type{KaplanMeier}, times, status)

Given a vector of times to events and a corresponding vector of indicators that
dictate whether each time is an observed event or is right censored, compute the
Kaplan-Meier estimate of the survivor function. Returns a [`KaplanMeier`](@ref)
object.
"""
function StatsBase.fit(::Type{KaplanMeier},
                       times::AbstractVector{T},
                       status::AbstractVector{<:Integer}) where {T}
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
    return _km(t, s)
end

function StatsBase.fit(::Type{KaplanMeier}, ev::EventTimeVector)
    isempty(ev) && throw(ArgumentError("the sample must be nonempty"))
    ev_sorted = sort(ev)
    return _km(ev_sorted.times, ev_sorted.status)
end

"""
    confint(km::KaplanMeier, α=0.05)

Compute the pointwise log-log transformed confidence intervals for the survivor
function as a vector of tuples.
"""
function StatsBase.confint(km::KaplanMeier, α::Float64=0.05)
    q = quantile(Normal(), 1 - α/2)
    return map(km.survival, km.stderr) do srv, se
        l = log(-log(srv))
        a = q * se / log(srv)
        exp(-exp(l - a)), exp(-exp(l + a))
    end
end
