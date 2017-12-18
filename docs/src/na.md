# Nelson-Aalen Estimator

The [Nelson-Aalen estimator](https://en.wikipedia.org/wiki/Kaplan-Meier_estimator)
is a nonparametric estimator of the cumulative hazard function, i.e. the probability of survival
beyond a given time.

The estimate is given by

```math
\hat{S}(t) = \sum_{i: t_i < t} \frac{d_i}{n_i}
```

where ``d_i`` is the number of observed events at time ``t_i`` and ``n_i`` is the
number of subjects at risk for the event just before time ``t_i``.

The pointwise standard error of the log of the survivor function can be computed
directly as the standard error or a Bernoulli random variable with `d_i` successes
from `n_i` samples:

```math
\text{SE}(\hat{S}(t)) = \sqrt{\sum_{i: t_i < t} \frac{d_i(n_i-d_i)}{n_i^3}}
```

## API

```@docs
Survival.NelsonAalen
StatsBase.fit(::Type{S},
              times::AbstractVector{T},
              status::AbstractVector{<:Integer}) where {S<:NonparametricEstimator, T}
StatsBase.confint(na::NelsonAalen, α::Float64=0.05)
```

## References

* Nelson, W. (1969). *Hazard plotting for incomplete failure data*.
  Journal of Quality Technology 1, 27–52.