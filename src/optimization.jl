function newton_raphson(fgh!,x; ρ = 0.5, c = 1e-4, tol = 1e-4, max_iter = 1000)
    grad= zeros(length(x))
    hes= zeros(length(x),length(x))
    y = 0.
    for i = 1:max_iter
        y = fgh!(x, grad, hes, true)
        search_dir = -(cholfact(Positive, hes)\grad)
        (norm(search_dir) > tol) || return x, y, grad, hes
        step_size = 1.
        while fgh!(x+search_dir*step_size, grad, hes, false) > y+c*step_size*(grad' * search_dir)
            step_size *= ρ
            step_size > 1e-10 || error("Linesearch failed! Problematic Hessian or Gradient?")
        end
        x = x + search_dir*step_size
    end
    error("Reached maximum number of iterations")
end
