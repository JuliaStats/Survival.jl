# preprocessed already:
# fs = index first deaths, ls = index last deaths,
# X is covariates, ξ is covariate covariate transpose
function _cox_h!(grad,hes, S, fs, ls, ξ , X, β, λ, alive, afterΘ, afterXΘ, afterξΘ, compute_derivatives)
    #compute relevant quantities for loglikelihood, score, fischer_info

    Xβ = X*β
    Θ = exp.(Xβ)
    after!(afterΘ,Θ)


    #compute loglikelihood, score, fischer_info
    y = 0.

    if compute_derivatives
        grad .= 0.
        hes .= 0.

        after!(afterXΘ, X.*Θ)
        after!(afterξΘ,ξ.*Θ)
        # preallocate
        Z = zeros(size(X,2))
        Ξ = zeros(size(X,2),size(X,2))
    end

    for i in 1:length(fs)
        for j in (fs[i]):(ls[i])
            ρ = (alive[j]-alive[fs[i]])/(alive[fs[i]]-alive[ls[i]+1])
            ϕ = afterΘ[fs[i]]-ρ*(afterΘ[fs[i]]-afterΘ[ls[i]+1])
            y -= Xβ[j] -log(ϕ)
            if compute_derivatives
                for k in eachindex(Z)
                    @inbounds Z[k] = afterXΘ[fs[i],k]-ρ*(afterXΘ[fs[i],k]-afterXΘ[ls[i]+1,k])
                end
                for k2 in 1:size(Ξ,2), k1 in 1:size(Ξ,1)
                    @inbounds Ξ[k1,k2] = afterξΘ[fs[i],k1,k2]-ρ*(afterξΘ[fs[i],k1,k2]-afterξΘ[ls[i]+1,k1,k2])
                end
                for k2 in 1:size(Ξ,2)
                    @inbounds grad[k2] -= X[j,k2]
                    @inbounds grad[k2] += Z[k2]/ϕ
                    for k1 in 1:size(Ξ,1)
                        @inbounds hes[k1, k2] += Ξ[k1,k2]/ϕ - Z[k1]*Z[k2]/ϕ^2
                    end
                end
            end
        end
    end
    y += λ*(β' * β)

    if compute_derivatives
        for k1 in 1:size(Ξ,1)
            grad[k1] +=  2*λ*β[k1]
            hes[k1,k1] +=  2*λ
        end
    end
    return y
end

function _coxph(X::AbstractArray{T}, S::AbstractVector; l2_cost = zero(T), kwargs...) where T
    ξ = zeros(size(X,1),size(X,2),size(X,2))
    for k2 in 1:size(ξ,3), k1 in 1:size(ξ,2), i in 1:size(ξ,1)
        @inbounds ξ[i,k1,k2] = X[i,k1]*X[i,k2]
    end

    # compute first and last!

    fs = find(firsts(S))
    ls = find(lasts(S))

    # do optimization

    alive = after(ones(Int64, length(S)))
    afterΘ = init_after(zeros(size(X,1)))
    afterXΘ = init_after(X.*zeros(size(X,1)))
    afterξΘ = init_after(ξ.*zeros(size(X,1)))

    f1 = (β, grad,hes) -> _cox_h!(grad,hes, S, fs, ls, ξ , X, β, l2_cost,alive, afterΘ, afterXΘ, afterξΘ, false)
    h1! = (β,grad,hes) -> _cox_h!(grad,hes, S, fs, ls, ξ , X, β, l2_cost,alive, afterΘ, afterXΘ, afterξΘ, true)
    β, neg_ll,grad, hes = newton_raphson(f1,h1!, zeros(size(X,2)); kwargs...)
    CoxModel(β, -neg_ll, -grad, hes)
end

function StatsBase.fit(::Type{CoxModel}, M::AbstractMatrix, y::AbstractVector; kwargs...)
    index_perm = sortperm(y)
    y_t = y[index_perm]
    M_t = M[index_perm, :]
    _coxph(M_t, y_t; kwargs...)
end

coxph(M, y; kwargs...) = fit(CoxModel, M, y; kwargs...)
