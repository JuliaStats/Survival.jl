# fs = index first deaths, ls = index last deaths,
# X is covariates, ξ is covariate covariate transpose

struct CoxAux{T}
    X::Array{T,2}
    ξ::Array{T,3}
    Xβ::Array{T,1}
    θ::Array{T,1}
    Xθ::Array{T,2}
    ξθ::Array{T,3}
    afterθ::Array{T,1}
    afterXθ::Array{T,2}
    afterξθ::Array{T,3}
    λ::T
    fs::Array{Int64,1}
    ls::Array{Int64,1}
    alive::Array{Int64,1}
end

function CoxAux(X::AbstractArray{T}, s::AbstractVector, l2_cost::T) where T
    ξ = zeros(T, size(X,1),size(X,2),size(X,2))
    for k2 in 1:size(ξ,3), k1 in 1:k2, i in 1:size(ξ,1)
        @inbounds ξ[i,k1,k2] = X[i,k1]*X[i,k2]
    end

    Xβ = zeros(T,size(X,1))
    θ = zeros(T,size(X,1))
    Xθ = zeros(T, size(X))
    ξθ =  zeros(T, size(ξ))
    afterθ = init_after(θ)
    afterXθ = init_after(Xθ)
    afterξθ = init_after(ξθ)

    fs = find(firsts(s))
    ls = find(lasts(s))
    alive = after(ones(Int64, length(s)))
    return CoxAux(X, ξ, Xβ, θ, Xθ, ξθ, afterθ, afterXθ, afterξθ, l2_cost, fs, ls, alive)
end

function update_cox!(c, β, compute_derivatives)
    A_mul_B!(c.Xβ, c.X, β)
    c.θ .= exp.(c.Xβ)
    after!(c.afterθ, c.θ)

    if compute_derivatives
        c.Xθ .= c.X .* c.θ
        c.ξθ .= c.ξ .* c.θ
        after!(c.afterXθ, c.Xθ)
        after!(c.afterξθ, c.ξθ)
    end
    return
end

function _cox_fgh!(β, grad,hes, c::CoxAux{T}, compute_derivatives)::T where T
    #get relevant quantities to compute negative loglikelihood, gradient and hessian

    update_cox!(c, β, compute_derivatives)

    X, ξ, Xβ, afterθ, afterXθ, afterξθ, λ, fs, ls, alive  =
        c.X, c.ξ, c.Xβ, c.afterθ, c.afterXθ, c.afterξθ, c.λ, c.fs, c.ls, c.alive


    #compute negative loglikelihood, gradient and hessian
    y = zero(T)

    if compute_derivatives
        grad .= zero(T)
        hes .= zero(T)

        # preallocate
        Z = zeros(T, length(β))
    end

    @inbounds for i in 1:length(fs)
        for j in (fs[i]):(ls[i])
            ρ = (alive[j]-alive[fs[i]])/(alive[fs[i]]-alive[ls[i]+1])
            ϕ = afterθ[fs[i]]-ρ*(afterθ[fs[i]]-afterθ[ls[i]+1])
            y -= Xβ[j] -log(ϕ)
            if compute_derivatives
                for k in eachindex(Z)
                    Z[k] = afterXθ[fs[i],k]-ρ*(afterXθ[fs[i],k]-afterXθ[ls[i]+1,k])
                end
                for k2 in 1:length(β)
                    grad[k2] += Z[k2]/ϕ - X[j,k2]
                    for k1 in 1:k2
                        Ξ = afterξθ[fs[i],k1,k2]-ρ*(afterξθ[fs[i],k1,k2]-afterξθ[ls[i]+1,k1,k2])
                        hes[k1, k2] += Ξ/ϕ - Z[k1]*Z[k2]/ϕ^2
                    end
                end
            end
        end
    end
    y += λ*(β' * β)

    if compute_derivatives
        for k1 in 1:length(β)
            grad[k1] +=  2*λ*β[k1]
            hes[k1,k1] +=  2*λ
        end
        for k2 in 1:length(β)
            for k1 in (k2+1):length(β)
                hes[k1, k2] = hes[k2, k1]
            end
        end
    end
    return y
end

function _coxph(X::AbstractArray{T}, s::AbstractVector; l2_cost = zero(T), kwargs...) where T
    c = CoxAux(X, s, l2_cost)

    fgh! = (β,grad,hes, compute_derivatives) ->
        _cox_fgh!(β, grad, hes, c, compute_derivatives)
    β, neg_ll,grad, hes = newton_raphson(fgh!, zeros(T, size(X,2)); kwargs...)
    CoxModel(β, -neg_ll, -grad, hes)
end

function StatsBase.fit(::Type{CoxModel}, M::AbstractMatrix, y::AbstractVector; kwargs...)
    index_perm = sortperm(y)
    X = M[index_perm, :]
    s = y[index_perm]
    _coxph(X, s; kwargs...)
end

coxph(M, y; kwargs...) = fit(CoxModel, M, y; kwargs...)
