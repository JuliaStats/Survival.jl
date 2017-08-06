function newton_raphson(f_,h_!,x; ρ = 0.5, c = 1e-4, tol = 1e-4, max_iter = 1000)
    grad= zeros(length(x))
    hes= zeros(length(x),length(x))
    y = 0.
    for i = 1:max_iter
        y = h_!(x, grad, hes)
        search_dir = -(cholfact(Positive, hes)\grad)
        (norm(search_dir) > tol) || break
        step_size = 1.
        while f_(x+search_dir*step_size, grad, hes) > y+c*step_size*(grad' * search_dir)
            step_size *= ρ
            step_size > 1e-10 || error("Linesearch failed! Problematic Hessian or Gradient?")
        end
        x = x + search_dir*step_size
    end
    return x, y, grad, hes
end
