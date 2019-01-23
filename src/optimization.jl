PositiveFactorizations.floattype(::Type{Union{T,Missing}}) where {T} = float(T)

function newton_raphson(fgh!, x::AbstractArray{T}; ρ=one(T)/2, c=1e-4, tol=1e-4, max_iter=1000) where T
    n = length(x)
    R = Base.promote_op(/, T, T)
    grad = zeros(R, n)
    hes = zeros(R, n, n)
    y = zero(R)
    for i = 1:max_iter
        y = fgh!(x, grad, hes, true)
        search_dir = -(cholesky(Positive, hes) \ grad)
        norm(search_dir) > tol || return x, y, grad, hes
        step_size = one(T)
        while fgh!(x + search_dir * step_size, grad, hes, false) > y + c * step_size * grad'search_dir
            step_size *= ρ
            step_size > 1e-10 || error("Line search failed! Problematic Hessian or Gradient?")
        end
        x .+= search_dir * step_size
    end
    throw(ConvergenceException(max_iter))
end
