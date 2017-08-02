struct CoxModel{T <: Real} <: RegressionModel
    model::AbstractString
    β::Array{T,1}
    loglik::T
    score::Array{T,1}
    fischer_info::Array{T,2}
end

function StatsBase.coeftable(obj::CoxModel)
    β, hes = obj.β,obj.fischer_info
    #rownms = coefnames(obj.mf)[2:end]
    se = sqrt.(diag(pinv(hes)))
    z_score = β./se
    pvalues = 2*cdf(Normal(),-abs.(z_score))
    coefmat = CoefTable(hcat([β, se, z_score, pvalues]...),
    ["Estimate", "Std.Error", "z value", "Pr(>|z|)"], ["x$i" for i in 1:length(β)], 4)
end

StatsBase.coef(obj::CoxModel) = obj.β
