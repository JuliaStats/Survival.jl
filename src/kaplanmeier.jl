"""
    KaplanMeier{S,T}

An immutable type containing survivor function estimates computed
using the Kaplan-Meier method.
The type has the following fields:

* `events`: An [`EventTable`](@ref) summarizing the times and events
  used to compute the estimates. The time values are of type `T`.
* `survival`: Estimate of the survival probability at each time. Values
  are of type `S`.
* `stderr`: Standard error of the log survivor function at each time.
  Values are of type `S`.

Use `fit(KaplanMeier, ...)` to compute the estimates as `Float64`
values and construct this type.
Alternatively, `fit(KaplanMeier{S}, ...)` may be used to request a
particular value type `S` for the estimates.
"""
struct KaplanMeier{S,T} <: NonparametricEstimator
    events::EventTable{T}
    survival::Vector{S}
    stderr::Vector{S}
end

"""
To computes the survival estimation at time `t` of a KaplanMeier fit you can use.

```julia
km = fit(KaplanMeier, ...)
t = 5
km(t) # evaluates the estimator at time 5
```
"""
function (km::KaplanMeier)(t::Real)
    time_points = km.events.time
    survival = km.survival
    if t < time_points[1]
        return estimator_start(typeof(km))
    else
        id = findlast(x -> x <= t, time_points) 
        return survival[id]
    end
end

estimator_eltype(::Type{<:KaplanMeier{S}}) where {S} = S
estimator_eltype(::Type{KaplanMeier}) = Float64

estimator_start(T::Type{<:KaplanMeier}) = oneunit(estimator_eltype(T))
stderr_start(T::Type{<:KaplanMeier}) = zero(estimator_eltype(T))

estimator_update(::Type{<:KaplanMeier}, es, dᵢ, nᵢ) = es * (1 - dᵢ // nᵢ)
stderr_update(::Type{<:KaplanMeier}, gw, dᵢ, nᵢ) = gw + dᵢ // (nᵢ * (nᵢ - dᵢ))

"""
    confint(km::KaplanMeier; level=0.05)

Compute the pointwise log-log transformed confidence intervals for the survivor
function as a vector of tuples.
"""
function StatsAPI.confint(km::KaplanMeier; level::Real=0.05)
    q = quantile(Normal(), 1 - level/2)
    return map(km.survival, km.stderr) do srv, se
        l = log(-log(srv))
        a = q * se / log(srv)
        exp(-exp(l - a)), exp(-exp(l + a))
    end
end

"""
    fit(KaplanMeier, times, status) -> KaplanMeier

Given a vector of times to events and a corresponding vector of indicators that
denote whether each time is an observed event or is right censored, compute the
Kaplan-Meier estimate of the survivor function.
"""
StatsAPI.fit(::Type{KaplanMeier}, times, status)

"""
    fit(KaplanMeier, ets) -> KaplanMeier

Compute the Kaplan-Meier estimate of the survivor function from a vector of
[`EventTime`](@ref) values.
"""
StatsAPI.fit(::Type{KaplanMeier}, ets)
