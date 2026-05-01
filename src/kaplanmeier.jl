"""
    KaplanMeier{S,T}

Immutable type containing survivor function estimates computed using the Kaplan-Meier
method. The estimates are evaluated at each unique observed time and are stepwise
constant between them. The fields are parallel vectors with one entry per unique
observed time:

* `events::EventTable{T}`: the [`EventTable`](@ref) summarizing the times and events
  used to compute the estimates.
* `survival::Vector{S}`: estimate of the survivor function at each time.
* `stderr::Vector{S}`: Greenwood's-formula standard error of the log of the survivor
  function at each time. Note that this is the standard error of `log(survival)`, not
  of `survival` itself; this form is convenient for log-log–transformed confidence
  intervals (see [`confint`](@ref)).

Use `fit(KaplanMeier, ...)` to compute the estimates as `Float64` values and construct
this type. Alternatively, `fit(KaplanMeier{S}, ...)` may be used to request a particular
value type `S` for the estimates (e.g. `Float32`).
"""
struct KaplanMeier{S,T} <: NonparametricEstimator
    events::EventTable{T}
    survival::Vector{S}
    stderr::Vector{S}
end

estimator_eltype(::Type{<:KaplanMeier{S}}) where {S} = S
estimator_eltype(::Type{KaplanMeier}) = Float64

estimator_start(T::Type{<:KaplanMeier}) = oneunit(estimator_eltype(T))
stderr_start(T::Type{<:KaplanMeier}) = zero(estimator_eltype(T))

estimator_update(::Type{<:KaplanMeier}, es, dᵢ, nᵢ) = es * (1 - dᵢ // nᵢ)
stderr_update(::Type{<:KaplanMeier}, gw, dᵢ, nᵢ) = gw + dᵢ // (nᵢ * (nᵢ - dᵢ))

"""
    confint(km::KaplanMeier; level::Real=0.95) -> Vector{Tuple{Float64,Float64}}

Compute pointwise confidence intervals for the survivor function at each time in
`km.events.time` and return them as a vector of `(lower, upper)` tuples. `level` is the
nominal coverage of the two-sided interval (default 0.95).

The intervals are constructed on the log-log scale, i.e. by forming a symmetric Wald
interval for `log(-log(S(t)))` and back-transforming. This guarantees that the bounds
lie in `[0, 1]` and tends to give better small-sample coverage than the plain Wald
interval on `S(t)` itself.
"""
function StatsAPI.confint(km::KaplanMeier; level::Real=0.95)
    q = quantile(Normal(), (1 + level)/2)
    return map(km.survival, km.stderr) do srv, se
        l = log(-log(srv))
        a = q * se / log(srv)
        exp(-exp(l - a)), exp(-exp(l + a))
    end
end

"""
    fit(::Type{KaplanMeier}, times::AbstractVector, status::AbstractVector) -> KaplanMeier
    fit(::Type{KaplanMeier}, ets::AbstractVector{<:EventTime}) -> KaplanMeier
    fit(::Type{KaplanMeier}, et::EventTable) -> KaplanMeier

Compute the Kaplan-Meier estimate of the survivor function. The input may be supplied
as parallel vectors of times and event-status indicators, as a vector of
[`EventTime`](@ref) values, or as a precomputed [`EventTable`](@ref). When `times` and
`status` are passed separately they must have the same length; `status` values are
converted to `Bool` (true ⇒ observed event, false ⇒ right censoring).

`fit(KaplanMeier{S}, ...)` may be used to request a particular value type `S` for the
estimates, e.g. `Float32`. The default is `Float64`.
"""
StatsAPI.fit(::Type{KaplanMeier}, times, status)
