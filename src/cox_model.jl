struct CoxModel{T <: Real} <: RegressionModel
    aux::CoxAux{T}
    β::Array{T,1}
    loglik::T
    score::Array{T,1}
    fischer_info::Array{T,2}
    vcov::Array{T,2}
end

function StatsBase.coeftable(obj::CoxModel)
    β = coef(obj)
    se = stderr(obj)
    z_score = β./se
    pvalues = 2*cdf(Normal(), -abs.(z_score))
    coefmat = CoefTable(hcat([β, se, z_score, pvalues]...),
    ["Estimate", "Std.Error", "z value", "Pr(>|z|)"], ["x$i" for i in 1:length(β)], 4)
end

function Base.show(io::IO, obj::CoxModel)
    println(io, "$(typeof(obj)):\n\nCoefficients:\n", coeftable(obj))
end

StatsBase.coef(obj::CoxModel) = obj.β

StatsBase.loglikelihood(obj::CoxModel) = obj.loglik
StatsBase.nullloglikelihood(obj::CoxModel{T}) where T = -_cox_f(obj.β*zero(T), obj.aux)

StatsBase.nobs(obj::CoxModel) = size(obj.aux.X, 1)

StatsBase.dof(obj::CoxModel) = length(obj.β)

StatsBase.vcov(obj::CoxModel) = obj.vcov

StatsBase.stderr(obj::CoxModel) = sqrt.(diag(vcov(obj)))
