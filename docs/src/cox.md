# Cox Proportional Hazards Model

The [Cox proportional hazards model](https://en.wikipedia.org/wiki/Proportional_hazards_model)
is a semiparametric regression model used to fit survival models without knowing the
distribution. It is based on the assumption that covariates affect the hazard function
multiplicatively. That is,

```math
\lambda(t | X_i) = \lambda_0(t) \exp(X_i \cdot \beta)
```

where ``\lambda(t|X_i)`` is the estimated hazard for sample ``i``, ``\lambda_0`` is the
baseline hazard, ``X_i`` is the vector of covariates for sample ``i``, and ``\beta`` is
the vector of coefficients in the model.

## API

```@docs
StatsAPI.fit(::Type{CoxModel}, M::AbstractMatrix, y::AbstractVector; kwargs...)
```

## References

* Cox, D. R. (1972). *Regression models and life tables (with discussion)*.
  Journal of the Royal Statistical Society, Series B, 34:187–220.
