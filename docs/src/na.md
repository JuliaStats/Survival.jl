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
```

## References

* Kaplan, E. L., and Meier, P. (1958). *Nonparametric Estimation from Incomplete
  Observations*. Journal of the American Statistical Association, 53(282), 457-481.
  doi:10.2307/2281868

* Greenwood, M. (1926). *A Report on the Natural Duration of Cancer*. Reports on
  Public Health and Medical Subjects. London: Her Majesty's Stationery Office.
  33, 1-26.
