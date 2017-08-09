function newton_raphson(fgh!, x::AbstractArray{T}; ρ = one(T)/2, c = 1e-4, tol = 1e-4, max_iter = 1000) where T
    grad = zeros(T, length(x))
    hes = zeros(T, length(x),length(x))
    y = zero(T)
    for i = 1:max_iter
        y = fgh!(x, grad, hes, true)
        search_dir = -(cholfact(Positive, hes)\grad)
        (norm(search_dir) > tol) || return x, y, grad, hes
        step_size = one(T)
        while fgh!(x+search_dir*step_size, grad, hes, false) > y+c*step_size*(grad' * search_dir)
            step_size *= ρ
            step_size > 1e-10 || error("Linesearch failed! Problematic Hessian or Gradient?")
        end
        x .+= search_dir*step_size
    end
    error("Reached maximum number of iterations")
end
