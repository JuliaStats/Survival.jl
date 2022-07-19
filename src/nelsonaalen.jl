"""
    NelsonAalen

An immutable type containing cumulative hazard function estimates computed
using the Nelson-Aalen method.
The type has the following fields:

* `times`: Distinct event times
* `nevents`: Number of observed events at each time
* `ncensor`: Number of right censored events at each time
* `natrisk`: Size of the risk set at each time
* `chaz`: Estimate of the cumulative hazard at each time
* `stderr`: Standard error of the cumulative hazard
Use `fit(NelsonAalen, ...)` to compute the estimates and construct
this type.
"""
struct NelsonAalen{T<:Real} <: NonparametricEstimator
    times::Vector{T}
    nevents::Vector{Int}
    ncensor::Vector{Int}
    natrisk::Vector{Int}
    chaz::Vector{Float64}
    stderr::Vector{Float64}
end

estimator_start(::Type{NelsonAalen}) = 0.0  # Estimator starting point
stderr_start(::Type{NelsonAalen}) = 0.0 # StdErr starting point

estimator_update(::Type{NelsonAalen}, es, dᵢ, nᵢ) = es + dᵢ / nᵢ # Estimator update rule
stderr_update(::Type{NelsonAalen}, gw, dᵢ, nᵢ) = gw + dᵢ * (nᵢ - dᵢ) / (nᵢ^3) # StdErr update rule

"""
    confint(na::NelsonAalen, α=0.05)

Compute the pointwise confidence intervals for the cumulative hazard
function as a vector of tuples.
"""
function StatsAPI.confint(na::NelsonAalen, α::Float64=0.05)
    q = quantile(Normal(), 1 - α/2)
    return map(na.chaz, na.stderr) do srv, se
        srv - q * se, srv + q * se
    end
end

"""
    fit(NelsonAalen, times, status) -> NelsonAalen

Given a vector of times to events and a corresponding vector of indicators that
denote whether each time is an observed event or is right censored, compute the
Nelson-Aalen estimate of the cumulative hazard rate function.
"""
StatsAPI.fit(::Type{NelsonAalen}, times, status)

"""
    fit(NelsonAalen, ets) -> NelsonAalen

Compute the Nelson-Aalen estimate of the cumulative hazard rate function from a
vector of [`EventTime`](@ref) values.
"""
StatsAPI.fit(::Type{NelsonAalen}, ets)
