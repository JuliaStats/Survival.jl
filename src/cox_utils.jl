# Compute first and last
firsts(s) = [s[t].status && (t==1 || s[t] > s[t-1]) for t = 1:length(s)]
lasts(s) = [s[t].status && (t==length(s) || s[t+1] > s[t]) for t = 1:length(s)]


# Compute "complementary cumulative sum" on first dimension
function after(v)
    afterv = init_after(v)
    after!(afterv,v)
    return afterv
end


function init_after(v::AbstractArray{T}) where T
    newsize = collect(size(v))
    newsize[1] += 1
    afterv = zeros(T,newsize...)
    return afterv
end

function after!{N,T}(afterv, v::AbstractArray{T,N})
    trails = fill(:,N-1)
    cumsum!(@view(afterv[(end-1):-1:1,trails...]),@view(v[end:-1:1,trails...]),1)
    return
end

struct CoxSum{T,N}
    values::Array{T,N}
    sums::Array{T,N}
    chunks::Array{T,N}
    tails::Array{T,N}
end

function CoxSum(values::Array{T,N}, fs, ls) where {T, N}
    sums = init_after(values)
    chunks = zeros(T, length(fs), Base.tail(size(values))...)
    tails = zeros(T, length(fs), Base.tail(size(values))...)
    CoxSum(values, sums, chunks, tails)
end

function update!(cs::CoxSum, fs, ls)
    after!(cs.sums, cs.values)
    for I in CartesianRange(Base.tail(indices(cs.sums)))
        for i in 1:length(fs)
            cs.chunks[i, I] = cs.sums[fs[i],I]-cs.sums[ls[i]+1,I]
            cs.tails[i, I] = cs.sums[fs[i], I]
        end
    end
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
    alive::Array{Int64,1}
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

    alive = after(ones(Int64, length(s)))
    return CoxAux(X, ξ, Xβ, θ, Xθ, ξθ,  l2_cost, fs, ls, alive)
end

function update_cox!(c, β, compute_derivatives)
    A_mul_B!(c.Xβ, c.X, β)
    c.θ.values .= exp.(c.Xβ)
    update!(c.θ, c.fs, c.ls)

    if compute_derivatives
        c.Xθ.values .= c.X .* c.θ.values
        c.ξθ.values .= c.ξ .* c.θ.values
        update!(c.Xθ, c.fs, c.ls)
        update!(c.ξθ, c.fs, c.ls)
    end
    return
end


# function get_chunks(M::AbstractArray{T}, fs, ls) where T
#     Mnew = zeros(T, length(fs), Base.tail(size(M))...)
#     for I in CartesianRange(Base.tail(indices(M)))
#         s = zero(T)
#         for i in size(M,1):-1:1
#             Mnew[i, I] = Mnew[i+1, I] + M[i,I]
#         end
#     end
#     return Mnew
# end
