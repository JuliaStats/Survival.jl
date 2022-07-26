# Event Times

A crucial concept in survival analysis is the time elapsed between some landmark and
a particular event of interest.
As an example, say you're running a clinical trial to investigate the efficacy of a
new anticonvulsant.
You may be interested in the time from the start of therapy to the first epileptic
seizure for each patient.
But if a patient dies or otherwise goes off study before they have a seizure, you'll
assume that a seizure would have occurred eventually, but you don't know when exactly.
In this case the event time is *right censored*; the true event time is unknown, all
you know is that it exceeds the observed time.

A dedicated type is provided to conveniently store right censored data.

```@docs
Survival.EventTime
```

## Summarizing Event Times

Given times to an event of interest and indications of whether the observations are
right censored, we can construct a table of the unique times along with the number of
events of interest, the number of censored observations, and the size of the risk set
at each time.
This information is used for computing other estimates, e.g. of the survivor and
cumulative hazard functions.

```@docs
Survival.EventTable
```
