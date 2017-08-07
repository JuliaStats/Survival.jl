# Compute first and last
firsts(s) = [s[t].status && (t==1 || s[t] > s[t-1]) for t = 1:length(s)]
lasts(s) = [s[t].status && (t==length(s) || s[t+1] > s[t]) for t = 1:length(s)]


struct CoxSum{T,N}
    values::Array{T,N}
    sums::Array{T,N}
    chunks::Array{T,N}
    tails::Array{T,N}
end

function CoxSum(values::Array{T,N}, fs, ls) where {T, N}
    sums = zeros(T, size(values,1)+1, Base.tail(size(values))...)
    chunks = zeros(T, length(fs), Base.tail(size(values))...)
    tails = zeros(T, length(fs), Base.tail(size(values))...)
    CoxSum(values, sums, chunks, tails)
end


# fs = index first deaths, ls = index last deaths,
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

function CoxAux(X::AbstractArray{T}, s::AbstractVector, l2_cost::T) where T
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

    return CoxAux(X, ξ, Xβ, θ, Xθ, ξθ,  l2_cost, fs, ls)
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
