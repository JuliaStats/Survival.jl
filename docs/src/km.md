# Kaplan-Meier Estimator

The [Kaplan-Meier estimator](https://en.wikipedia.org/wiki/Kaplan-Meier_estimator)
is a nonparametric estimator of the survivor function, i.e. the probability of survival
beyond a given time.

The estimate is given by

```math
\hat{S}(t) = \prod_{i: t_i < t} \left( 1 - \frac{d_i}{n_i} \right)
```

where ``d_i`` is the number of observed events at time ``t_i`` and ``n_i`` is the
number of subjects at risk for the event just before time ``t_i``.

The pointwise standard error of the log of the survivor function can be computed
using Greenwood's formula:

```math
\text{SE}(\log \hat{S}(t)) = \sqrt{\sum_{i: t_i < t} \frac{d_i}{n_i (n_i - d_i)}}
```

## Plotting

Survival.jl provides an [AlgebraOfGraphics.jl](https://github.com/MakieOrg/AlgebraOfGraphics.jl)
extension for creating Kaplan-Meier plots. The extension is loaded automatically when both
Survival and AlgebraOfGraphics are available.

```@example km
using Survival, AlgebraOfGraphics, CairoMakie, RDatasets

lung = dataset("survival", "lung")
lung.Event = lung.Status .== 2
lung.SexLabel = replace.(string.(lung.Sex), "1" => "Male", "2" => "Female")
first(lung, 5)
```

### Basic Kaplan-Meier Plot

Map time and event indicator columns to the first two positional arguments of
`kaplanmeier()`:

```@example km
plt = data(lung) * mapping("Time", "Event") * (kaplanmeier() + censorticks())
draw(plt)
```

### Grouped Kaplan-Meier Plot

Use the `color` mapping to stratify by a grouping variable.
The `marker` mapping is applied only to `censorticks()` so that the censor
marks get group-specific marker shapes:

```@example km
plt = data(lung) * mapping("Time", "Event", color="SexLabel") *
    (kaplanmeier() + mapping(marker="SexLabel") * censorticks())
draw(plt)
```

### Adding a Risk Table

Call `add_risktable!` on the result of `draw()` to add a number-at-risk table
below the plot:

```@example km
plt = data(lung) * mapping("Time", "Event", color="SexLabel") *
    (kaplanmeier() + mapping(marker="SexLabel") * censorticks())
fg = draw(plt)
add_risktable!(fg)
fg
```

### Faceted by Column

Use the `col` mapping to display groups in separate panels.
`add_risktable!` works with faceted layouts as well:

```@example km
plt = data(lung) * mapping("Time", "Event", col="SexLabel") *
    (kaplanmeier() + censorticks())
fg = draw(plt)
add_risktable!(fg)
fg
```

### Custom Confidence Level

Use the `level` keyword to change the confidence level (default 0.95):

```@example km
plt = data(lung) * mapping("Time", "Event") * (kaplanmeier(level=0.9) + censorticks())
draw(plt)
```

### Without Confidence Interval

Pass `interval=nothing` to show only the survival curve:

```@example km
plt = data(lung) * mapping("Time", "Event") * kaplanmeier(interval=nothing)
draw(plt)
```

## API

```@docs
Survival.KaplanMeier
StatsAPI.fit(::Type{KaplanMeier}, ::Any, ::Any)
StatsAPI.fit(::Type{KaplanMeier}, ::Any)
StatsAPI.confint(::KaplanMeier)
Survival.kaplanmeier
Survival.censorticks
Survival.add_risktable!
```

## References

* Kaplan, E. L., and Meier, P. (1958). *Nonparametric Estimation from Incomplete
  Observations*. Journal of the American Statistical Association, 53(282), 457-481.
  doi:10.2307/2281868

* Greenwood, M. (1926). *A Report on the Natural Duration of Cancer*. Reports on
  Public Health and Medical Subjects. London: Her Majesty's Stationery Office.
  33, 1-26.
