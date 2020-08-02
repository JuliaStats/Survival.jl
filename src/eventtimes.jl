#-- EventTime ---------------





## Type constructors


"""
    EventTime{T}

Immutable object containing the real-valued time to an event as well as an indicator of
whether the time corresponds to an observed event (`true`) or right censoring (`false`).
"""
struct EventTime{T<:Real} <: AbstractEventTime
    time::T
    status::Bool
end

EventTime(time::T) where {T<:Real} = EventTime{T}(time, true)

## Overloaded Base functions

Base.eltype(::EventTime{T}) where {T} = T
Base.show(io::IO, ev::EventTime) = print(io, ev.time, ifelse(ev.status, "", "+"))

Base.convert(T::Type{<:Real}, ev::EventTime) = convert(T, ev.time)
Base.convert(T::Type{EventTime}, x::Real) = EventTime(x)

function Base.isless(t1::EventTime, t2::EventTime)
    if t1.time == t2.time
        # When two EventTimes have the same observed time, we compare the event
        # status. Observed events compare less than censored events, since if the
        # censored event were to occur then it would happen at or after the given
        # time (by definition).
        return isevent(t1) && iscensored(t2)
    else
        return isless(t1.time, t2.time)
    end
end

## New functions

iscensored(ev::EventTime) = !ev.status
isevent(ev::EventTime) = ev.status


## StatsModels compatibility

StatsModels.concrete_term(t::Term, xs::AbstractVector{<:EventTime}, ::Nothing) =
    StatsModels.ContinuousTerm(t.sym, first(xs), first(xs), first(xs), first(xs))
Base.copy(et::EventTime) = et















############## COMPETING EVENTS

## Type constructors

"""
`CompetingEventTime{T,S}(time, status, eventofinterest, censoringevent)`

Immutable object containing the real-valued time to an event as well as a status variable
indicating whether the time corresponds to the event of interest, a competing event, or censoring.
It also stores the values of `status` that correspond to the event of interest and censored observations.

By default, `status` is an integer valued variable with `1` being the event of interest and `0` being censored values.
When using the default values, it is possible to supply only `time` and `status`.
When status can take on other types like `String` or `Symbol`, the full constructor is required.

```
CompetingEventTime(1.7,1,1,0)
CompetingEventTime(1.7,1)   # same as above
CompetingEventTime(1.7,1; eventofinterest=5, censoringevent=-1)   # Int keyword args supported
CompetingEventTime(1.2,"death","relapse","censor")  # full constructor
CompetingEventTime(1.2,:death,:relapse,:censor) 
```
"""
struct CompetingEventTime{T<:Real,S}
    time::T
    status::S
    eventofinterest::S
    censoringevent::S
end

# CompetingEventTime(time::T) where {T<:Real} = CompetingEventTime{T,Int64}(time, 1, 1, 0)
# helper function for probably the more common case (status being Int)
# for all other cases, require user to use full 4 arg constructor
CompetingEventTime(time::T,status::S=1; eventofinterest::S=one(S), censoringevent::S=zero(S)) where {T<:Real,S<:Int} = CompetingEventTime{T,S}(time, status, eventofinterest, censoringevent)


## New functions
iscensored(ev::CompetingEventTime) = ev.status == ev.censoringevent
iseventofinterest(ev::CompetingEventTime) = ev.status == ev.eventofinterest
isevent(ev::CompetingEventTime) = !iscensored(ev)
iscompetingevent(ev::CompetingEventTime) = !iscensored(ev) && !iseventofinterest(ev)

eventtype(::CompetingEventTime{T,S}) where {T,S} = S

## Overloaded Base functions

Base.eltype(::CompetingEventTime{T,S}) where {T,S} = T
function Base.show(io::IO, ev::CompetingEventTime{T,S}) where {T,S} 
    # print(io, "CompetingEventTime{",T,",",S,"}(",ev.time, ifelse(iseventofinterest(ev), "", "+"),ifelse(iscensored(ev),"",string(" (",ev.status,")")),")")
    printstring = string(ev.time, ifelse(iseventofinterest(ev), "", "+"))
    if iscompetingevent(ev)
        printstring = string(printstring," (",ev.status,")")
    end
    print(io, printstring)
end


Base.convert(T::Type{<:Real}, ev::CompetingEventTime) = convert(T, ev.time)
Base.convert(T::Type{CompetingEventTime}, x::Real) = CompetingEventTime(x)

function Base.isless(t1::CompetingEventTime, t2::CompetingEventTime)
    if t1.time == t2.time
        # When two EventTimes have the same observed time, we compare the event
        # status. Observed events compare less than censored events, since if the
        # censored event were to occur then it would happen at or after the given
        # time (by definition). Two tied events are considered to have no proper order
        # (all comparisons return false).
        return isevent(t1) && iscensored(t2)
    else
        return isless(t1.time, t2.time)
    end
end



swapeventofinterest(ev::CompetingEventTime{T,S},neweoi::S) where {T,S} = CompetingEventTime{T,S}(ev.time,ev.status,neweoi,ev.censoringevent)
swapcensoringevent(ev::CompetingEventTime{T,S},newcensor::S) where {T,S} = CompetingEventTime{T,S}(ev.time,ev.status,ev.eventofinterest,newcensor)
Base.copy(ev::CompetingEventTime) = ev

## StatsModels compatibility
# StatsModels.concrete_term(t::Term, xs::AbstractVector{<:CompetingEventTime}, ::Nothing) =
#     StatsModels.ContinuousTerm(t.sym, first(xs), first(xs), first(xs), first(xs))
