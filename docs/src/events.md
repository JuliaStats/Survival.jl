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

To conveniently store right censored data, two types are provided for convenience.

```@docs
Survival.EventTime
Survival.EventTimeVector
```
