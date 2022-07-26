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

EventTime(time, status=true) = EventTime{typeof(time)}(time, Bool(status))

## Overloaded Base functions

Base.eltype(::Type{EventTime{T}}) where {T} = T
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

## StatsModels compatibility

StatsModels.concrete_term(t::Term, xs::AbstractVector{<:EventTime}, ::Nothing) =
    StatsModels.ContinuousTerm(t.sym, first(xs), first(xs), first(xs), first(xs))

Base.copy(et::EventTime) = et

#####
##### `EventTable`
#####

"""
    EventTable{T}

Immutable object summarizing the unique observed event times, including the number of
events, the number of censored observations, and the number remaining at risk for each
unique time.

This type implements the Tables.jl interface for tables, which means that `EventTable`
objects can be seamlessly converted to other tabular types such as `DataFrame`s.

    EventTable(eventtimes)

Construct an `EventTable` from an array of [`EventTime`](@ref) values.

    EventTable(time, status)

Construct an `EventTable` from an array of time values and an array of event status
indicators.
"""
struct EventTable{T}
    time::Vector{T}
    nevents::Vector{Int}
    ncensored::Vector{Int}
    natrisk::Vector{Int}
end

function EventTable(ets)
    T = eltype(eltype(ets))
    isempty(ets) && return EventTable{T}(T[], Int[], Int[], Int[])
    ets = issorted(ets) ? ets : sort(ets)  # re-binding, input is unaffected
    _droptimezero!(ets)
    return _eventtable(ets)
end

function EventTable(time, status)
    ntimes = length(time)
    nstatus = length(status)
    if ntimes != nstatus
        throw(DimensionMismatch("number of event statuses does not match number of " *
                                "event times; got $nstatus and $ntimes, respectively"))
    end
    T = eltype(time)
    ntimes == 0 && return EventTable{T}(T[], Int[], Int[], Int[])
    ets = map(EventTime, time, status)
    issorted(ets) || sort!(ets)
    _droptimezero!(ets)
    return _eventtable(ets)
end

function _droptimezero!(ets)
    # Assumptions about the input:
    #   - iterates `EventTime`s
    #   - sorted ascending by elements' `.time` fields
    i = findfirst(et -> !iszero(et.time), ets)
    start = firstindex(ets)
    if i !== nothing && i > start
        deleteat!(ets, start:(start + i - 1))
    end
    return ets
end

function _eventtable(ets)
    # Assumptions about the input:
    #   - nonempty
    #   - time 0 is not included
    T = typeof(first(ets).time)
    outlen = _nuniquetimes(ets)

    nobs = length(ets)
    dᵢ::Int = 0                   # Number of observed events at time t
    cᵢ::Int = 0                   # Number of censored events at time t
    nᵢ::Int = nobs                # Number remaining at risk at time t

    times = Vector{T}(undef, outlen)            # The set of unique event times
    nevents = Vector{Int}(undef, outlen)        # Total observed events at each time
    ncensor = Vector{Int}(undef, outlen)        # Total censored events at each time
    natrisk = Vector{Int}(undef, outlen)        # Number at risk at each time

    t_prev = zero(T)
    outind = 1

    @inbounds begin
        for et in ets
            t = et.time
            s = et.status
            # Aggregate over tied times
            if t == t_prev
                dᵢ += s
                cᵢ += !s
                continue
            elseif !iszero(t_prev)
                times[outind] = t_prev
                nevents[outind] = dᵢ
                ncensor[outind] = cᵢ
                natrisk[outind] = nᵢ
                outind += 1
            end
            nᵢ -= dᵢ + cᵢ
            dᵢ = s
            cᵢ = !s
            t_prev = t
        end

        # We need to do this one more time to capture the last time
        # since everything in the loop is lagged
        times[outind] = t_prev
        nevents[outind] = dᵢ
        ncensor[outind] = cᵢ
        natrisk[outind] = nᵢ
    end

    return EventTable{eltype(times)}(times, nevents, ncensor, natrisk)
end

function _nuniquetimes(ets)
    # Assumptions about the input:
    #   - nonempty
    #   - iterates `EventTime`s
    #   - sorted ascending by elements' `.time` fields
    t_prev = first(ets).time
    n = 1
    for et in Iterators.drop(ets, 1)
        t = et.time
        if t != t_prev
            n += 1
            t_prev = t
        end
    end
    return n
end

Base.copy(et::EventTable{T}) where {T} =
    EventTable{T}(copy(et.time), copy(et.nevents), copy(et.ncensored), copy(et.natrisk))

Base.:(==)(a::EventTable, b::EventTable) =
    a.time == b.time && a.nevents == b.nevents &&
    a.ncensored == b.ncensored && a.natrisk == b.natrisk

# Tables.jl integration

Tables.istable(::Type{<:EventTable}) = true

Tables.columnaccess(::Type{<:EventTable}) = true

_rowtype(T::Type{<:EventTable}) =
    NamedTuple{fieldnames(T),Tuple{map(eltype, fieldtypes(T))...}}

_rowtype(et::EventTable) = _rowtype(typeof(et))

Tables.schema(et::EventTable) = Tables.Schema(_rowtype(et))

function Tables.rows(et::EventTable)
    NT = _rowtype(et)
    nr = length(et.time)
    nc = fieldcount(NT)
    return (@inbounds(NT(ntuple(i -> getfield(et, i)[j], nc))) for j in 1:nr)
end
