#####
##### `fit` for non-parametric estimators
#####

function StatsAPI.fit(::Type{E}, et::EventTable) where {E<:NonparametricEstimator}
    outlen = length(et.time)
    outlen == 0 && throw(ArgumentError("can't fit `$(nameof(E))` on 0 observations"))
    T = eltype(et.time)
    S = estimator_eltype(E)
    estimator = Vector{S}(undef, outlen)
    stderror = Vector{S}(undef, outlen)
    es = estimator_start(E)
    gw = stderr_start(E)
    @inbounds for i in 1:outlen
        dᵢ = et.nevents[i]
        nᵢ = et.natrisk[i]
        es = estimator_update(E, es, dᵢ, nᵢ)
        gw = stderr_update(E, gw, dᵢ, nᵢ)
        estimator[i] = es
        stderror[i] = sqrt(gw)
    end
    return Core.apply_type(Base.typename(E).wrapper, S, T)(et, estimator, stderror)
end

StatsAPI.fit(::Type{E}, args...) where {E<:NonparametricEstimator} =
    fit(E, EventTable(args...))
