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


If your data contain competing risks, you should use a `CompetingEventTime`. This type allows you to specify which type of event occurred and what your event of interest is.

```@docs
Survival.CompetingEventTime
```

## API

```@docs
Survival.iscensored
Survival.isevent
Survival.iseventofinterest
Survival.eventtype
Survival.iscompetingevent(::CompetingEventTime)
Survival.swapeventofinterest(::CompetingEventTime)
Survival.swapcensoringevent(::CompetingEventTime)
```