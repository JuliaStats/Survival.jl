function _cox_f(β, c::CoxAux{T})::T where T
    grad = zeros(T, length(β))
    hes = zeros(T, length(β), length(β))
    _cox_fgh!(β, grad, hes, c::CoxAux{T}, false)
end

function _cox_fgh!(β, grad, hes, c::CoxAux{T}, compute_derivatives)::T where T
    #get relevant quantities to compute negative loglikelihood, gradient and hessian

    update_cox!(c, β, compute_derivatives)

    X, ξ, Xβ, θ, Xθ, ξθ, λ, fs, ls  =
        c.X, c.ξ, c.Xβ, c.θ, c.Xθ, c.ξθ, c.λ, c.fs, c.ls


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
            ρ = (j-fs[i])/(ls[i] + 1 - fs[i])
            ϕ = θ.tails[i]-ρ * θ.chunks[i]
            y -= Xβ[j] - log(ϕ)
            if compute_derivatives
                for k in eachindex(Z)
                    Z[k] = Xθ.tails[i, k] - ρ * Xθ.chunks[i, k]
                end
                for k2 in 1:length(β)
                    grad[k2] += Z[k2]/ϕ - X[j, k2]
                    for k1 in 1:k2
                        Ξ = ξθ.tails[i, k1, k2] - ρ * ξθ.chunks[i, k1, k2]
                        hes[k1, k2] += Ξ/ϕ - Z[k1] * Z[k2]/ϕ^2
                    end
                end
            end
        end
    end
    y += λ*(β' * β)

    if compute_derivatives
        for k1 in 1:length(β)
            grad[k1] += 2*λ*β[k1]
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

_coxph(X::AbstractArray{T}, s::AbstractVector; kwargs...) where {T<:Integer} =
    _coxph(Float64.(X), s::AbstractVector; kwargs...)

function _coxph(X::AbstractArray{T}, s::AbstractVector; l2_cost = zero(T), kwargs...) where T
    c = CoxAux(X, s, l2_cost)

    fgh! = (β, grad, hes, compute_derivatives) ->
        _cox_fgh!(β, grad, hes, c, compute_derivatives)
    β, neg_ll, grad, hes = newton_raphson(fgh!, zeros(T, size(X,2)); kwargs...)
    CoxModel(c, β, -neg_ll, -grad, hes, pinv(hes))
end

"""
    fit(::Type{CoxModel}, M::AbstractMatrix, y::AbstractVector; kwargs...)

Given a matrix M of predictors and a corresponding vector of events, compute the
Cox proportional hazard model estimate of coefficients. Returns a [`CoxModel`](@ref)
object.
"""
function StatsBase.fit(::Type{CoxModel}, M::AbstractMatrix, y::AbstractVector; kwargs...)
    index_perm = sortperm(y)
    X = M[index_perm, :]
    s = y[index_perm]
    _coxph(X, s; kwargs...)
end

coxph(M, y; kwargs...) = fit(CoxModel, M, y; kwargs...)
