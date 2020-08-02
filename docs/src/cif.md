# Cumulative Incidence Functions

In the presence of competing risks, the Kaplan-Meier survival function may not be the appropriate estimate of survival probability. The cumulative incidence function (CIF) can be used to study the probability of failure by a particular type of event instead. 

In the presence of competing risks, `1-KM`, where `KM` is the Kaplan-Meier estimate,
is uninterpretable and is a biased estimate of the failure probability. The
cumulative incidence estimator of Kalbfleisch and Prentice (1980) is a function
of the hazards of both the event of interest and the competing event, and provides
an unbiased estimate of the failure probability.

The estimator of the cumulative incidence for event ``k`` is given by
```math
\hat{I}_k(t) = \sum_{i: t_i < t} \hat{S}(t_{i-1}) \frac{d_{k,i}}{n_i}
```
where ``d_{k,i}`` are the events of interest ``k`` at time ``t_i``,
``n_i`` are the individuals at risk at that time,
and ``\\hat{S}(t_{i-1})`` is the usual Kaplan-Meier estimate of survival. Standard errors are computed using the Delta method.

This package implements a nonparametric `CumulativeIncidence` estimator. The best way to use the `CumulativeIncidence` type is to cast your data into a vector of `CompetingEventTime` objects and then pass it to fit `fit(CumulativeIncidence,...)`.


## API

```@docs
Survival.CumulativeIncidence
StatsBase.fit(::Type{CumulativeIncidence},tte::AbstractVector{CompetingEventTime{T,S}}) where {T<:Real,S}
StatsBase.confint(::CumulativeIncidence, ::Float64)
```

## References

* Kalbfleisch, J. D., and Prentice, R. L. (1980). *The Statistical Analysis of
  Failure Time Data*. New York, NY: John Wiley.