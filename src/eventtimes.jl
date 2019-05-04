#-- EventTime ---------------

## Type constructors

"""
    EventTime{T}

Immutable object containing the real-valued time to an event as well as an indicator of
whether the time corresponds to an observed event (`true`) or right censoring (`false`).
"""
struct EventTime{T<:Real}
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

isevent(ev::EventTime) = ev.status
iscensored(ev::EventTime) = !ev.status


# Concrete term to play along StatsModels

struct EventTimeTerm <: AbstractTerm
    sym::Symbol
end

function concrete_term(t::Term, xs::AbstractVector{<:EventTime}, ::Nothing)
    EventTimeTerm(t.sym)
end

modelcols(t::EventTimeTerm, d::NamedTuple) = d[t.sym]