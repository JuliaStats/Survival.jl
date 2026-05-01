"""
    NelsonAalen{S,T}

Immutable type containing cumulative-hazard function estimates computed using the
Nelson-Aalen method. The estimates are evaluated at each unique observed time and are
stepwise constant between them. The fields are parallel vectors with one entry per
unique observed time:

* `events::EventTable{T}`: the [`EventTable`](@ref) summarizing the times and events
  used to compute the estimates.
* `chaz::Vector{S}`: estimate of the cumulative hazard at each time.
* `stderr::Vector{S}`: standard error of the cumulative hazard at each time, computed
  as the binomial-variance form ``\\sqrt{\\sum_i d_i(n_i - d_i)/n_i^3}``.

Use `fit(NelsonAalen, ...)` to compute the estimates as `Float64` values and construct
this type. Alternatively, `fit(NelsonAalen{S}, ...)` may be used to request a particular
value type `S` for the estimates (e.g. `Float32`).
"""
struct NelsonAalen{S,T} <: NonparametricEstimator
    events::EventTable{T}
    chaz::Vector{S}
    stderr::Vector{S}
end

estimator_eltype(::Type{<:NelsonAalen{S}}) where {S} = S
estimator_eltype(::Type{NelsonAalen}) = Float64

estimator_start(T::Type{<:NelsonAalen}) = zero(estimator_eltype(T))
stderr_start(T::Type{<:NelsonAalen}) = zero(estimator_eltype(T))

estimator_update(::Type{<:NelsonAalen}, es, dᵢ, nᵢ) = es + dᵢ // nᵢ
stderr_update(::Type{<:NelsonAalen}, gw, dᵢ, nᵢ) = gw + dᵢ * (nᵢ - dᵢ) // (nᵢ^3)

"""
    confint(na::NelsonAalen; level::Real=0.95) -> Vector{Tuple{Float64,Float64}}

Compute pointwise Wald confidence intervals for the cumulative hazard function at each
time in `na.events.time` and return them as a vector of `(lower, upper)` tuples.
`level` is the nominal coverage of the two-sided interval (default 0.95). The lower
bound is not clipped at zero, so it may be negative for small estimates.
"""
function StatsAPI.confint(na::NelsonAalen; level::Real=0.95)
    q = quantile(Normal(), (1 + level)/2)
    return map(na.chaz, na.stderr) do srv, se
        srv - q * se, srv + q * se
    end
end

"""
    fit(::Type{NelsonAalen}, times::AbstractVector, status::AbstractVector) -> NelsonAalen
    fit(::Type{NelsonAalen}, ets::AbstractVector{<:EventTime}) -> NelsonAalen
    fit(::Type{NelsonAalen}, et::EventTable) -> NelsonAalen

Compute the Nelson-Aalen estimate of the cumulative hazard function. The input may be
supplied as parallel vectors of times and event-status indicators, as a vector of
[`EventTime`](@ref) values, or as a precomputed [`EventTable`](@ref). When `times` and
`status` are passed separately they must have the same length; `status` values are
converted to `Bool` (true ⇒ observed event, false ⇒ right censoring).

`fit(NelsonAalen{S}, ...)` may be used to request a particular value type `S` for the
estimates, e.g. `Float32`. The default is `Float64`.
"""
StatsAPI.fit(::Type{NelsonAalen}, times, status)
