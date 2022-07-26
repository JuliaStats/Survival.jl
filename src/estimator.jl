#####
##### `fit` for non-parametric estimators
#####

function StatsAPI.fit(::Type{S}, et::EventTable) where {S<:NonparametricEstimator}
    outlen = length(et.time)
    outlen == 0 && throw(ArgumentError("can't fit `$(nameof(S))` on 0 observations"))
    estimator = Vector{Float64}(undef, outlen)
    stderror = Vector{Float64}(undef, outlen)
    es = estimator_start(S)
    gw = stderr_start(S)
    @inbounds for i in 1:outlen
        dᵢ = et.nevents[i]
        nᵢ = et.natrisk[i]
        es = estimator_update(S, es, dᵢ, nᵢ)
        gw = stderr_update(S, gw, dᵢ, nᵢ)
        estimator[i] = es
        stderror[i] = sqrt(gw)
    end
    return S{eltype(et.time)}(et, estimator, stderror)
end

StatsAPI.fit(::Type{S}, args...) where {S<:NonparametricEstimator} =
    fit(S, EventTable(args...))
