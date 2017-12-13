# Some utility functions to compute Cox regression

# Compute first and last
firsts(s) = (s[t].status && (t==1 || s[t] > s[t-1]) for t = 1:length(s))
lasts(s) = (s[t].status && (t==length(s) || s[t+1] > s[t]) for t = 1:length(s))


struct CoxSum{T,N}
    values::Array{T,N}
    sums::Array{T,N}
    chunks::Array{T,N}
    tails::Array{T,N}
end

function CoxSum(values::Array{T,N}, fs, ls) where {T, N}
    tailsize = Base.tail(size(values))
    sums = zeros(T, size(values,1)+1, tailsize...)
    chunks = zeros(T, length(fs), tailsize...)
    CoxSum(values, sums, chunks, copy(chunks))
end


# fs, ls = indices first and last events in a simultaneaous group
# if there are no ties, fs == ls
# X is covariates, ξ is covariate covariate transpose

struct CoxAux{T}
    X::Array{T,2}
    ξ::Array{T,3}
    Xβ::Array{T,1}
    θ::CoxSum{T,1}
    Xθ::CoxSum{T,2}
    ξθ::CoxSum{T,3}
    λ::T
    fs::Array{Int64,1}
    ls::Array{Int64,1}
end

function CoxAux(X::AbstractArray{T}, s::AbstractVector, l2_cost) where T
    fs = find(firsts(s))
    ls = find(lasts(s))
    ξ = zeros(T, size(X,1),size(X,2),size(X,2))
    for k2 in 1:size(ξ,3), k1 in 1:k2, i in 1:size(ξ,1)
        @inbounds ξ[i,k1,k2] = X[i,k1]*X[i,k2]
    end

    Xβ = zeros(T,size(X,1))
    θ = CoxSum(zeros(T,size(X,1)), fs, ls)
    Xθ = CoxSum(zeros(T,size(X)), fs, ls)
    ξθ =  CoxSum(zeros(T, size(ξ)), fs, ls)

    return CoxAux(X, ξ, Xβ, θ, Xθ, ξθ,  T(l2_cost), fs, ls)
end

function update_cox!(c::CoxAux, β, compute_derivatives)
    A_mul_B!(c.Xβ, c.X, β)
    c.θ.values .= exp.(c.Xβ)
    update_cox!(c.θ, c.fs, c.ls)

    if compute_derivatives
        c.Xθ.values .= c.X .* c.θ.values
        c.ξθ.values .= c.ξ .* c.θ.values
        update_cox!(c.Xθ, c.fs, c.ls)
        update_cox!(c.ξθ, c.fs, c.ls)
    end
    return
end

function update_cox!(cs::CoxSum, fs, ls)
    @inbounds for I in CartesianRange(Base.tail(indices(cs.values)))
        if (length(I) < 2) || (I[1] <= I[2])
            for i in size(cs.values)[1]:-1:1
                cs.sums[i, I] = cs.sums[i+1, I] + cs.values[i, I]
            end
            for i in 1:length(fs)
                cs.chunks[i, I] = cs.sums[fs[i],I]-cs.sums[ls[i]+1,I]
                cs.tails[i, I] = cs.sums[fs[i], I]
            end
        end
    end
end

# Structure of Cox regression output

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
    pvalues = 2 .* cdf.(Normal(), -abs.(z_score))
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

#compute negative loglikelihood

function _cox_f(β, c::CoxAux{T})::T where T
    grad = zeros(T, length(β))
    hes = zeros(T, length(β), length(β))
    _cox_fgh!(β, grad, hes, c::CoxAux{T}, false)
end

#compute negative loglikelihood, gradient and hessian

function _cox_fgh!(β, grad, hes, c::CoxAux{T}, compute_derivatives)::T where T
    update_cox!(c, β, compute_derivatives)

    X, ξ, Xβ, θ, Xθ, ξθ, λ, fs, ls  =
        c.X, c.ξ, c.Xβ, c.θ, c.Xθ, c.ξθ, c.λ, c.fs, c.ls

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

StatsModels.supports_intercept(::Type{CoxModel}) = false

"""
    fit(::Type{CoxModel}, M::AbstractMatrix, y::AbstractVector; kwargs...)

Given a matrix `M` of predictors and a corresponding vector of events, compute the
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
