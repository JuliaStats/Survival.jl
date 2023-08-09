# Nelson-Aalen Estimator

The [Nelson-Aalen estimator](https://en.wikipedia.org/wiki/Nelson%E2%80%93Aalen_estimator)
is a nonparametric estimator of the cumulative hazard function.

The estimate is given by

```math
\hat{H}(t) = \sum_{i: t_i < t} \frac{d_i}{n_i}
```

where ``d_i`` is the number of observed events at time ``t_i`` and ``n_i`` is the
number of subjects at risk for the event just before time ``t_i``.

The pointwise standard error of the log of the survivor function can be computed
directly as the standard error or a Bernoulli random variable with ``d_i`` successes
from ``n_i`` samples:

```math
\text{SE}(\hat{H}(t)) = \sqrt{\sum_{i: t_i < t} \frac{d_i(n_i-d_i)}{n_i^3}}
```

## API

```@docs
Survival.NelsonAalen
Survival.NelsonAalen(t::Real)
StatsAPI.fit(::Type{NelsonAalen}, ::Any, ::Any)
StatsAPI.confint(::NelsonAalen)
```

## References

* Nelson, W. (1969). *Hazard plotting for incomplete failure data*.
  Journal of Quality Technology 1, 27–52.
