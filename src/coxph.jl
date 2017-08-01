function cox_f(S, fs, ls, ξ , X, β, λ, alive, afterΘ) # preprocessed already
    Xβ = X*β
    Θ = exp.(Xβ)
    after!(afterΘ, Θ)

    #compute loglikelihood
    y = 0.
    for i in 1:length(fs)
        for j in (fs[i]):(ls[i])
            ρ = (alive[j]-alive[fs[i]])/(alive[fs[i]]-alive[ls[i]+1])
            ϕ = afterΘ[fs[i]]-ρ*(afterΘ[fs[i]]-afterΘ[ls[i]+1])
            y -= Xβ[j] -log(ϕ)
        end
    end
    y += λ*(β' * β)
    return y
end

# preprocessed already:
# fs = index first deaths, ls = index last deaths,
# X is covariates, ξ is covariate covariate transpose
function cox_h!(grad,hes, S, fs, ls, ξ , X, β, λ, alive, afterΘ, afterXΘ, afterξΘ)
    #compute relevant quantities for loglikelihood, score, fischer_info

    Xβ = X*β
    Θ = exp.(Xβ)
    after!(afterΘ,Θ)
    after!(afterXΘ, X.*Θ)
    after!(afterξΘ,ξ.*Θ)

    #compute loglikelihood, score, fischer_info
    #From v0.6 remember to use . notation and do it inplace, possibly using views!
    y = 0.
    grad[:] = 0.
    hes[:] = 0.

    # preallocate
    Z = zeros(size(X,2))
    Ξ = zeros(size(X,2),size(X,2))

    for i in 1:length(fs)
        for j in (fs[i]):(ls[i])
            ρ = (alive[j]-alive[fs[i]])/(alive[fs[i]]-alive[ls[i]+1])
            ϕ = afterΘ[fs[i]]-ρ*(afterΘ[fs[i]]-afterΘ[ls[i]+1])
            for k in eachindex(Z)
                @inbounds Z[k] = afterXΘ[fs[i],k]-ρ*(afterXΘ[fs[i],k]-afterXΘ[ls[i]+1,k])
            end
            for k2 in 1:size(Ξ,2)
                for k1 in 1:size(Ξ,1)
                    @inbounds Ξ[k1,k2] = afterξΘ[fs[i],k1,k2]-ρ*(afterξΘ[fs[i],k1,k2]-afterξΘ[ls[i]+1,k1,k2])
                end
            end
            y -= Xβ[j] -log(ϕ)
            for k2 in 1:size(Ξ,2)
                @inbounds grad[k2] -= X[j,k2]
                @inbounds grad[k2] += Z[k2]/ϕ
                for k1 in 1:size(Ξ,1)
                    @inbounds hes[k1, k2] += Ξ[k1,k2]/ϕ - Z[k1]*Z[k2]/ϕ^2
                end
            end
        end
    end
    y += λ*(β' * β)

    for k1 in 1:size(Ξ,1)
        grad[k1] +=  2*λ*β[k1]
        hes[k1,k1] +=  2*λ
    end
    return y
end

function coxph(S::AbstractVector,X::AbstractArray; l2_cost = 0., kwargs...)
    ξ = zeros(size(X,1),size(X,2),size(X,2))
    for k2 in 1:size(ξ,3)
        for k1 in 1:size(ξ,2)
            for i in 1:size(ξ,1)
                @inbounds ξ[i,k1,k2] = X[i,k1]*X[i,k2]
            end
        end
    end

    # compute first and last!

    fs = find(firsts(S))
    ls = find(lasts(S))

    # do optimization

    alive = after(ones(Int64, length(S)))
    afterΘ = init_after(zeros(size(X,1)))
    afterXΘ = init_after(X.*zeros(size(X,1)))
    afterξΘ = init_after(ξ.*zeros(size(X,1)))

    f1 = (β) -> cox_f(S, fs, ls, ξ , X, β, l2_cost, alive, afterΘ)
    h1! = (β,grad,hes) -> cox_h!(grad,hes, S, fs, ls, ξ , X, β, l2_cost,alive, afterΘ, afterXΘ, afterξΘ)
    return newton_raphson(f1,h1!, zeros(size(X,2)); kwargs...)
end

function coxph(formula::Formula, data::DataFrame; l2_cost = 0., kwargs...)
    sorted_data = sort(data, cols = formula.lhs)
    M = DataFrames.ModelFrame(formula,sorted_data)
    S = convert(Array, M.df[:,1])
    model_matrix = DataFrames.ModelMatrix(M)
    X = model_matrix.m[:,2:size(model_matrix.m,2)]
    β, neg_ll,grad, hes =  coxph(S, X; l2_cost = l2_cost, kwargs...)
    rownms = coefnames(M)[2:end]
    se = sqrt.(diag(pinv(hes)))
    z_score = β./se
    pvalues = 2*cdf(Normal(),-abs.(z_score))
    coefmat = CoefTable(hcat([β, se, z_score, pvalues]...),
    ["Estimate", "Std.Error", "z value", "Pr(>|z|)"], rownms, 4)
    CoxModel("Cox; ", formula, coefmat, M, -neg_ll, -grad, hes)
end
