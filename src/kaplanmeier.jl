"""
    KaplanMeier

An immutable type containing survivor function estimates computed
using the Kaplan-Meier method.
The type has the following fields:

* `times`: Distinct event times
* `nevents`: Number of observed events at each time
* `ncensor`: Number of right censored events at each time
* `natrisk`: Size of the risk set at each time
* `survival`: Estimate of the survival probability at each time
* `stderr`: Standard error of the log survivor function at each time

Use `fit(KaplanMeier, ...)` to compute the estimates and construct
this type.
"""
struct KaplanMeier{T<:Real} <: NonparametricEstimator
    times::Vector{T}
    nevents::Vector{Int}
    ncensor::Vector{Int}
    natrisk::Vector{Int}
    survival::Vector{Float64}
    stderr::Vector{Float64}
end

estimator_start(::Type{KaplanMeier}) = 1.0  # Estimator starting point
stderr_start(::Type{KaplanMeier}) = 0.0 # StdErr starting point

estimator_update(::Type{KaplanMeier}, es, dᵢ, nᵢ) = es * (1 - dᵢ / nᵢ) # Estimator update rule
stderr_update(::Type{KaplanMeier}, gw, dᵢ, nᵢ) = gw + dᵢ / (nᵢ * (nᵢ - dᵢ)) # StdErr update rule

"""
    confint(km::KaplanMeier, α=0.05)

Compute the pointwise log-log transformed confidence intervals for the survivor
function as a vector of tuples.
"""
function StatsBase.confint(km::KaplanMeier, α::Float64=0.05)
    q = quantile(Normal(), 1 - α/2)
    return map(km.survival, km.stderr) do srv, se
        l = log(-log(srv))
        a = q * se / log(srv)
        exp(-exp(l - a)), exp(-exp(l + a))
    end
end

"""
    fit(KaplanMeier, times, status) -> KaplanMeier

Given a vector of times to events and a corresponding vector of indicators that
dictate whether each time is an observed event or is right censored, compute the
Kaplan-Meier estimate of the survivor function.
"""
StatsBase.fit(::Type{KaplanMeier}, times, status)

# Calculate the number at risk in each group and total number of
# events at the given time point 't', which may not be an event time
# in some of the KaplanMeier fits.
function update_km!(km::NTuple{N,KaplanMeier{T}}, t::T, natrisk::Vector,
                    nevents::Vector) where {N,T}

    natrisk .= 0
    nevents .= 0
    for (j,k) in enumerate(km)
        l = searchsortedfirst(k.times, t)
        if l > length(k.times)
            l -= 1
        else
            natrisk[j] = k.natrisk[l]
        end
        if k.times[l] == t
            nevents[j] = k.nevents[l]
        end
    end
end

function get_stat(od::Vector, va::Matrix, trend=nothing)
	p = length(od) - 1
    stat, dof = if isnothing(trend)
    	# Standard test
    	try
    		od[1:p]' * (va[1:p, 1:p] \ od[1:p]), p
		catch SingularException
			@warn("Singular covariance encountered in log rank test, using pseudo-inverse")
    		od[1:p]' * pinv(va[1:p, 1:p]) * od[1:p], p
		end
	else
		# Trend test
		dot(trend, od)^2 / (trend' * va * trend), 1
	end
	return stat, dof
end

"""
    logrank_test(km) -> NamedTuple

Conduct a logrank test for the null hypothesis that n population
survival functions are equal.  The parameters here are an arbitrary
number of fitted survival functions of type 'KaplanMeier'.
"""
function logrank_test(km::KaplanMeier...; method=:LogRank, psf=nothing, trend=nothing, fh=[1., 1.])

    if method == :FH && isnothing(psf)
        error("When method is 'FH', the pooled survival function 'psf' must be provided")
    end

    # We are comparing p groups
    p = length(km)

    # All distinct times across all groups
    tu = sort(unique(vcat([k.times for k in km]...)))

    # Total events per group
    obs = zeros(p)

    # Workspace
    natrisk, nevents, expval = zeros(p), zeros(p), zeros(p)

    # The variance/covariance matrix of the estimates
    va = zeros(p, p)

	# Accumulate the weights
    wtsum = 0.0
    wtcount = 0.0

    # Dictionary to locate the entry in psf that corresponds
    # to each unique time.
    psf_ix = Dict{Int,Int}()
    if !isnothing(psf)
        for i in eachindex(psf.times)
            psf_ix[psf.times[i]] = i
        end
    end

    for (i,t) in enumerate(tu)

        update_km!(km, t, natrisk, nevents)
        r = sum(natrisk)
        d = sum(nevents)

        wt = if method == :LogRank
                1.0
            elseif method == :WBG
                r
            elseif method == :TW
                sqrt(r)
            elseif method == :FH
                i = psf_ix[t]
                f = i > 1 ? psf.survival[i-1] : 1
                f^fh[1] * (1 - f)^fh[2]
            else
                error("Unknown weighting method in logrank_test")
            end
        wtsum += wt
        wtcount += 1

        # Update the expected values
        expval .+= d * wt .* natrisk ./ sum(natrisk)

        # Update the number of events
        obs .+= wt .* nevents

        # Update the variance/covariance matrix
        if r > 1
            f = d*(r-d) / (r^2*(r-1))
            for i1 in 1:p
                va[i1, i1] += wt^2 * f * natrisk[i1] * (r - natrisk[i1])
                for i2 in 1:i1-1
                    u = wt^2 * f * natrisk[i1] * natrisk[i2]
                    va[i1, i2] -= u
                    va[i2, i1] -= u
                end
            end
        end
    end

    # Normalize the weights
    scale = wtcount ./ wtsum
    obs .*= scale
    expval .*= scale
    va .*= scale.^2

    # Compute the test statistic
    od = obs - expval
    stat, dof = get_stat(od, va, trend)
    pvalue = 1 - cdf(Chisq(dof), stat)

    return (stat=stat, dof=dof, pvalue=pvalue, observed=obs, expected=expval, var=va)
end

"""
    logrank_test(km) -> NamedTuple

Test the null hypothesis that n population survival functions are
equal using stratified data.  The parameter km is a vector of vectors
of fitted survival functions of type 'KaplanMeier'.  The outer level
of this vector of vectors corresponds to distinct strata, i.e. km[1]
is the first stratum, km[2] is the second stratum, etc.  The inner
level corresponds to groups, i.e. km[2][3] is the third group in the
second stratum.  The null hypothesis is that the population survival
functions are equal across all groups within each stratum.
"""
function logrank_test(km::Vector{Vector{KaplanMeier}}; method=:LogRank,
                      psf=fill(nothing, length(km)), trend=nothing, fh=[1.0, 1.0])

    if method == :FH && (nothing in psf)
        error("If method is FH, the pooled survival function psf must be provided")
    end

    # Conduct a logrank test for each stratum
    st = [logrank_test(x...; method=method, psf=ps, fh=fh) for (x, ps) in zip(km, psf)]
    observed = sum([x.observed for x in st])
    expected = sum([x.expected for x in st])
    variance = sum([x.var for x in st])

    # Compute the test statistic
    od = observed - expected
    stat, dof = get_stat(od, variance, trend)
    pvalue = 1 - cdf(Chisq(dof), stat)

    return (stat=stat, dof=dof, pvalue=pvalue, observed=observed, expected=expected, var=variance)
end

# Return the indices where a new consecutive run of identical values
# occurs in the sorted array 'v'.  The first such index is always 1,
# and the final index length(v) + 1 is included as a sentinel.
function index_splits(v)
    jj = [i for i in 2:length(v) if v[i] != v[i-1]]
    prepend!(jj, 1)
    push!(jj, length(v)+1)
    return jj
end

# Split the 'time' and 'status' vectors according to the distinct
# values in the vector 'group', then fit a Kaplan Meier estimate of
# the survival function for each subset of data thus obtained.  Also
# returns the labels of the distinct groups, sorted compatibly with
# the returned vector of KaplanMeier values.
function make_km(time, status, group)

    ii = sortperm(group)
    time = time[ii]
    status = status[ii]
    group = group[ii]

    jj = index_splits(group)
    km = KaplanMeier[]
    for j in 1:length(jj)-1
        i1, i2 = jj[j], jj[j+1]
        k = fit(KaplanMeier, time[i1:i2-1], status[i1:i2-1])
        push!(km, k)
    end

    return km, unique(group)
end

# Return a vector of vectors of KaplanMeier values, with the outer
# index corresponding to strata and the inner index corresponding to
# groups.  Also returns the unique strata labels, and the unique group
# labels within each stratum.
function make_km(time, status, group, strata)

    ii = sortperm(strata)
    time = time[ii]
    status = status[ii]
    group = group[ii]

    strat = sort(unique(strata))
    groups = []
    jj = index_splits(strata)
    km = Vector{Vector{KaplanMeier}}()
    for j in 1:length(jj)-1
        i1, i2 = jj[j], jj[j+1]
        kx, g = make_km(time[i1:i2-1], status[i1:i2-1], group[i1:i2-1])
        push!(km, kx)
        push!(groups, g)
    end

    return km, strat, groups
end

"""
    logrank_test(times, status, group) -> NamedTuple

Test the null hypothesis that n population survival functions are
equal.  The parameters of the function are vectors containing time
(numeric), status (1=event, 0=censored), and group labels, with the
group labels defining the samples to be compared.
"""
function logrank_test(time, status, group; method=:LogRank, trend=nothing, fh=[1., 1.])

    if !(length(time) == length(status) == length(group))
        @warn("'time', 'status', and 'group' must have equal lengths")
    end

	# Get the survival function pooling over groups
    psf = if method == :FH
            fit(KaplanMeier, time, status)
        else
            nothing
        end

    km, g = make_km(time, status, group)

	trendvec = isnothing(trend) ? nothing : [trend[x] for x in g]

    return logrank_test(km...; method=method, psf=psf, trend=trendvec, fh=fh)
end

# Return a vector containing fitted Kaplan Meier values for the
# distinct values of strata (sorted).
function pooled_km(time::AbstractVector, status::AbstractVector, strata::AbstractVector)
    ii = sortperm(strata)
    strata = strata[ii]
    time = time[ii]
    status = status[ii]
    jj = index_splits(strata)
    ps = KaplanMeier[]
    for j in 1:length(jj)-1
        j1, j2 = jj[j], jj[j+1]-1
        k = fit(KaplanMeier, time[j1:j2], status[j1:j2])
        push!(ps, k)
    end
    return ps
end

"""
    logrank_test(times, status, group, strata; method, fh) -> NamedTuple

Test the null hypothesis that n population survival functions are
equal, using stratified data.  The parameters of the function are
vectors containing time (numeric), status (1=event, 0=censored), group
labels, and strata labels.  The null hypothesis being tested is that
the population survival functions are identical across groups within
each stratum.

The 'method' argument determines the weighting used in the test.
Options are :LogRank (all weight are 1, the default), :WBG
(Wilcoxon-Breslow-Gehan, weighting by the number at risk), :TW
(Tarone-Ware, weighting by the square root of the number at risk), and
:FH (Fleming-Harrington, weighting by a function of the pooled
survival function).  If method is :FH, then the 'fh' argument is used
to specify the exponents in the weight, which is S^fh[1] * (1 -
S)^fh[2], with S being the estimated survival function at the previous
time point (S is calculated by pooling over groups but calculating
separately by stratum).
"""
function logrank_test(time, status, group, strata; method=:LogRank, trend=nothing, fh=[1.0, 1.0])

    if !(length(time) == length(status) == length(group) == length(strata))
        @warn("'time', 'status', 'group', and 'strata' must have equal lengths")
    end

    km, _, g = make_km(time, status, group, strata)

    # Get the pooled survival function for each
    # stratum (pooling over groups within strata).
    psf = if method == :FH
            pooled_km(time, status, strata)
        else
            fill(nothing, length(km))
        end

	trendvec = nothing
	if !isnothing(trend)
		for j in eachindex(g)
			if !all(g[j] .== g[1])
				@warn("Some strata do not have observations from all groups")
			end
		end
		trendvec = [trend[x] for x in g[1]]
	end

    return logrank_test(km; method=method, psf=psf, trend=trendvec, fh=fh)
end
