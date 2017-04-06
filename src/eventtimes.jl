#-- EventTime ---------------

## Type constructors

"""
    EventTime{T}

Immutable object containing the real-valued time to an event as well as an indicator of
whether the time corresponds to an observed event or right censoring.
"""
struct EventTime{T<:Real}
    time::T
    status::Bool
end

EventTime(time::T) where {T<:Real} = EventTime{T}(time, false)

## Overloaded Base functions

Base.eltype(::EventTime{T}) where {T} = T
Base.show(io::IO, ev::EventTime) = print(io, ev.time, ifelse(ev.status, "", "+"))

Base.convert(T::Type{<:Real}, ev::EventTime) = convert(T, ev.time)
Base.convert(T::Type{EventTime}, x::Real) = EventTime(x)

## New functions

isevent(ev::EventTime) = ev.status
iscensored(ev::EventTime) = !ev.status


#-- EventTimeVector ---------

## Type constructors

"""
    EventTimeVector{T}

A vector containing event times and indicators of whether each time is an observed event or
right censoring.
"""
mutable struct EventTimeVector{T<:Real} <: AbstractVector{T}
    times::Vector{T}
    status::BitArray{1}

    function EventTimeVector{T}(times::Vector{T}, status::BitVector) where {T<:Real}
        if length(times) != length(status)
            throw(DimensionMismatch("event times and statuses must have the same length"))
        end
        return new{T}(times, status)
    end
end

EventTimeVector(times::Vector, status::AbstractVector{Bool}) =
    EventTimeVector{eltype(times)}(times, BitArray(status))

function EventTimeVector(evtimes::AbstractVector{EventTime{T}}) where {T}
    n = length(evtimes)
    t = Vector{T}(n)
    s = BitVector(n)
    @inbounds for (i, evtime) in enumerate(evtimes)
        t[i] = evtime.time
        s[i] = evtime.status
    end
    return EventTimeVector{T}(t, s)
end

EventTimeVector(t::AbstractVector{T}) where {T<:Real} = EventTimeVector{T}(t, falses(t))

## Overloaded Base functions

for f in [:length, :endof, :size, :eachindex, :isempty]
    @eval Base.$f(ev::EventTimeVector) = $f(ev.times)
end

Base.indices(ev::EventTimeVector) = (Base.OneTo(length(ev.times)),)

Base.iteratorsize(::EventTimeVector) = Base.HasShape()
Base.iteratoreltype(::EventTimeVector) = Base.HasEltype()
Base.IndexStyle(::EventTimeVector) = IndexLinear()
Base.IndexStyle(::Type{EventTimeVector{T}}) where {T} = IndexLinear()

Base.start(ev::EventTimeVector) = 1
Base.next(ev::EventTimeVector, state::Integer) = (ev[state], state + 1)
Base.done(ev::EventTimeVector, state::Integer) = state > length(ev)

Base.eltype(::EventTimeVector{T}) where {T} = EventTime{T}

Base.getindex(ev::EventTimeVector, i::Integer) = EventTime(ev.times[i], ev.status[i])
Base.getindex(ev::EventTimeVector, r::AbstractArray{<:Integer}) =
    EventTimeVector(ev.times[r], ev.status[r])
Base.getindex(ev::EventTimeVector, ::Colon) = copy(ev)

function Base.setindex!(ev::EventTimeVector, t::EventTime, i::Integer)
    @boundscheck checkbounds(ev.times, i)
    @inbounds begin
        ev.times[i] = t.time
        ev.status[i] = t.status
    end
    return t
end

function Base.setindex!(ev::EventTimeVector, t::EventTime, r::AbstractArray{<:Integer})
    @boundscheck checkbounds(ev.times, r)
    @inbounds begin
        ev.times[r] = t.time
        ev.status[r] = t.status
    end
    return t
end

Base.broadcast(::typeof(isevent), ev::EventTimeVector) = ev.status
Base.broadcast(::typeof(iscensored), ev::EventTimeVector) = .!ev.status

Base.any(::typeof(isevent), ev::EventTimeVector) = any(ev.status)
Base.all(::typeof(isevent), ev::EventTimeVector) = all(ev.status)
Base.any(::typeof(iscensored), ev::EventTimeVector) = any(!, ev.status)
Base.all(::typeof(iscensored), ev::EventTimeVector) = all(!, ev.status)

function Base.copy(ev::EventTimeVector{T}) where {T}
    n = length(ev)
    t = copy!(Array{T}(n), ev.times)
    s = copy!(BitVector(n), ev.status)
    return EventTimeVector{T}(t, s)
end

Base.sort(ev::EventTimeVector) = ev[sortperm(ev.times)]

function Base.sort!(ev::EventTimeVector)
    p = sortperm(ev.times)
    permute!(ev.times, p)
    permute!(ev.status, p)
    return ev
end

function Base.push!(ev::EventTimeVector, t::EventTime)
    push!(ev.times, t.time)
    push!(ev.status, t.status)
    return ev
end

for f in [:append!, :prepend!]
    @eval function Base.$f(ev1::EventTimeVector, ev2::EventTimeVector)
        $f(ev1.times, ev2.times)
        $f(ev1.status, ev2.status)
        return ev1
    end
end

Base.vcat(ev::EventTimeVector, evs::EventTimeVector...) =
    EventTimeVector(vcat(ev.times, getfield.(evs, :times)...),
                    vcat(ev.status, getfield.(evs, :status)...))
