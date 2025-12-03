# Some utility functions to compute Cox regression

promote_nonmissing(::Type{T}) where {T} = typeof(one(T) / one(T))
promote_nonmissing(::Type{Union{T,Missing}}) where {T} = promote_nonmissing(T)

struct CoxSum{T,N}
    values::Array{T,N}
    sums::Array{T,N}
    chunks::Array{T,N}
    tails::Array{T,N}
end

function CoxSum(values::Array{T}, fs, ls) where T
    tailsize = Base.tail(size(values))
    sums = zeros(T, size(values, 1) + 1, tailsize...)
    chunks = zeros(T, length(fs), tailsize...)
    CoxSum(values, sums, chunks, copy(chunks))
end

# fs, ls = indices first and last events in a simultaneous group
# if there are no ties, fs == ls
# X is covariates, ξ is covariate covariate transpose

struct CoxAux{T}
    X::Matrix{T}
    ξ::Array{T,3}
    Xβ::Vector{T}
    θ::CoxSum{T,1}
    Xθ::CoxSum{T,2}
    ξθ::CoxSum{T,3}
    λ::T
    fs::Vector{Int64}
    ls::Vector{Int64}
end

function CoxAux(X::AbstractArray{T}, s::AbstractVector, l2_cost) where T
    R = promote_nonmissing(T)
    ξ = zeros(R, size(X, 1), size(X, 2), size(X, 2))
    @inbounds for k2 in 1:size(ξ, 3), k1 in 1:k2, i in 1:size(ξ,1)
        ξ[i,k1,k2] = X[i,k1]*X[i,k2]
    end
    fs = findall(t->s[t].status && (t == 1 || s[t] > s[t-1]), 1:length(s))
    ls = findall(t->s[t].status && (t == length(s) || s[t+1] > s[t]), 1:length(s))
    Xβ = zeros(R, size(X, 1))
    θ = CoxSum(zeros(R, size(X, 1)), fs, ls)
    Xθ = CoxSum(zeros(R, size(X)), fs, ls)
    ξθ = CoxSum(zeros(R, size(ξ)), fs, ls)
    return CoxAux{R}(X, ξ, Xβ, θ, Xθ, ξθ, R(l2_cost), fs, ls)
end

function update_cox!(c::CoxAux, β, compute_derivatives)
    mul!(c.Xβ, c.X, β)
    c.θ.values .= exp.(c.Xβ)
    update_cox!(c.θ, c.fs, c.ls)
    if compute_derivatives
        c.Xθ.values .= c.X .* c.θ.values
        c.ξθ.values .= c.ξ .* c.θ.values
        update_cox!(c.Xθ, c.fs, c.ls)
        update_cox!(c.ξθ, c.fs, c.ls)
    end
    nothing
end

function update_cox!(cs::CoxSum, fs, ls)
    @inbounds for I in CartesianIndices(Base.tail(axes(cs.values)))
        if length(I) < 2 || I[1] <= I[2]
            for i = size(cs.values, 1):-1:1
                cs.sums[i,I] = cs.sums[i+1,I] + cs.values[i,I]
            end
            for i = 1:length(fs)
                cs.chunks[i,I] = cs.sums[fs[i],I] - cs.sums[ls[i]+1,I]
                cs.tails[i,I] = cs.sums[fs[i],I]
            end
        end
    end
    nothing
end

# Structure of Cox regression output

struct CoxModel{T<:Real} <: RegressionModel
    aux::CoxAux{T}
    β::Vector{T}
    loglik::T
    score::Vector{T}
    fischer_info::Matrix{T}
    vcov::Matrix{T}
end

function StatsAPI.coeftable(obj::CoxModel)
    β = coef(obj)
    se = stderror(obj)
    z_score = β ./ se
    pvalues = 2 .* cdf.(Normal(), -abs.(z_score))
    CoefTable(hcat(β, se, z_score, pvalues),
              ["Estimate", "Std.Error", "z value", "Pr(>|z|)"],
              map(i->string("x", i), 1:length(β)), 4)
end

function Base.show(io::IO, model::CoxModel)
    ct = coeftable(model)
    println(io, typeof(model))
    println(io)
    println(io, "Coefficients:")
    show(io, ct)
end

StatsAPI.coef(obj::CoxModel) = obj.β

StatsAPI.modelmatrix(obj::CoxModel) = obj.aux.X

StatsAPI.loglikelihood(obj::CoxModel) = obj.loglik

StatsAPI.nullloglikelihood(obj::CoxModel) = -_cox_f(zero(coef(obj)), obj.aux)

StatsAPI.nobs(obj::CoxModel) = size(modelmatrix(obj), 1)

StatsAPI.dof(obj::CoxModel) = length(coef(obj))

StatsAPI.dof_residual(obj::CoxModel) = nobs(obj) - dof(obj)

StatsAPI.vcov(obj::CoxModel) = obj.vcov

StatsAPI.stderror(obj::CoxModel) = sqrt.(diag(vcov(obj)))

function StatsAPI.confint(obj::CoxModel; level::Real=0.95)
    β = coef(obj)
    se = stderror(obj)
    q = quantile(Normal(), (1 - level) / 2)
    ci_upr = β .- se .* q
    ci_lwr = β .+ se .* q
    return hcat(ci_lwr, ci_upr)
end

#compute negative loglikelihood

function _cox_f(β, c::CoxAux{T})::T where T
    grad = zeros(T, length(β))
    hes = zeros(T, length(β), length(β))
    _cox_fgh!(β, nothing, nothing, c)
end

#compute negative loglikelihood, gradient and hessian

function _cox_fgh!(β, grad, hes, c::CoxAux{T}) where T
    update_cox!(c, β, (hes !== nothing) | (grad !== nothing))
    @compat (; X, ξ, Xβ, θ, Xθ, ξθ, λ, fs, ls) = c
    y = zero(T) / one(T)

    if hes !== nothing
        fill!(hes, 0)
    end
    if grad !== nothing
        fill!(grad, 0)
    end
    if grad !== nothing || hes !== nothing
        # preallocate
        Z = zeros(typeof(y), length(β))
    end

    @inbounds for i in 1:length(fs)
        for j in fs[i]:ls[i]
            ρ = (j - fs[i]) / (ls[i] + 1 - fs[i])
            ϕ = θ.tails[i] - ρ * θ.chunks[i]
            y -= Xβ[j] - log(ϕ)
            if hes !== nothing || grad !== nothing
                for k in eachindex(Z)
                    Z[k] = Xθ.tails[i,k] - ρ * Xθ.chunks[i,k]
                end
                for k2 in 1:length(β)
                    if grad !== nothing
                        grad[k2] += Z[k2] / ϕ - X[j,k2]
                    end
                    if hes !== nothing
                        for k1 in 1:k2
                            Ξ = ξθ.tails[i,k1,k2] - ρ * ξθ.chunks[i,k1,k2]
                            hes[k1,k2] += Ξ / ϕ - Z[k1] * Z[k2] / ϕ^2
                        end
                    end
                end
            end
        end
    end
    y += λ*(β'β)
    if hes !== nothing || grad !== nothing
        for k1 in 1:length(β)
            if grad !== nothing
                grad[k1] += 2*λ*β[k1]
            end
            if hes !== nothing
                hes[k1,k1] += 2*λ
            end
        end
        if hes !== nothing
            for k2 in 1:length(β)
                for k1 in k2+1:length(β)
                    hes[k1,k2] = hes[k2,k1]
                end
            end
        end
    end
    return y
end

_coxph(X::AbstractArray{<:Integer}, s::AbstractVector; tol, l2_cost) = _coxph(float(X), s; tol=tol, l2_cost=l2_cost)

function _coxph(X::AbstractArray{T}, s::AbstractVector; l2_cost, tol) where T
    R = promote_nonmissing(T)
    c = CoxAux(X, s, l2_cost)
    β₀ = zeros(R, size(X, 2))
    fgh! = NLSolversBase.TwiceDifferentiable(NLSolversBase.only_fgh!((f, G, H, x)->_cox_fgh!(x, G, H, c)), β₀)
    optim_alg = Optim.NewtonTrustRegion(; delta_min = 0.0)
    optim_options = Optim.Options(; g_abstol = tol)
    optim_state = Optim.initial_state(optim_alg, optim_options, fgh!, β₀)
    res = Optim.optimize(fgh!, β₀, optim_alg, optim_options, optim_state)
    rescode = Optim.termination_code(res)
    @show rescode
    if rescode != Optim.TerminationCode.GradientNorm && rescode != Optim.TerminationCode.NoXChange && rescode != Optim.TerminationCode.NoObjectiveChange
        error(LazyString("Calculation of estimate of Cox proportional hazard model failed: Optimization stopped with termination code `", rescode, "`."))
    end
    β = Optim.minimizer(res)
    @assert β == optim_state.x
    neg_ll = Optim.minimum(res)
    @assert neg_ll == optim_state.f_x
    grad = optim_state.g_x
    hes = optim_state.H_x
    return CoxModel{R}(c, β, -neg_ll, -grad, hes, pinv(hes))
end

StatsModels.drop_intercept(::Type{CoxModel}) = true

"""
    fit(::Type{CoxModel}, M::AbstractMatrix, y::AbstractVector; kwargs...)

Given a matrix `M` of predictors and a corresponding vector of events, compute the
Cox proportional hazard model estimate of coefficients. Returns a `CoxModel`
object.
"""
function StatsAPI.fit(::Type{CoxModel}, M::AbstractMatrix, y::AbstractVector; tol=1e-4, l2_cost=0)
    index_perm = sortperm(y)
    X = M[index_perm,:]
    s = y[index_perm]
    _coxph(X, s; tol=tol, l2_cost=l2_cost)
end

coxph(M, y; kwargs...) = fit(CoxModel, M, y; kwargs...)
