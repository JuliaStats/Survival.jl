module SurvivalAlgebraOfGraphicsExt

using Survival
using Survival: KaplanMeier, EventTable, fit, confint
using AlgebraOfGraphics
using AlgebraOfGraphics: transformation, ProcessedLayer, ProcessedLayers, dictionary, verbatim
using AlgebraOfGraphics.Makie: Makie, Stairs, Band, Scatter, Axis, to_value, text!, hidexdecorations!, linkxaxes!, ylims!

import Survival: kaplanmeier, censorticks, add_risktable!

# ── Helper: stepped band coordinates ────────────────────────────────────────

# Duplicate points to create stepped band coordinates for use with `Band`,
# which otherwise linearly interpolates between points.
function step_band_coords(t, lo, hi)
    n = length(t)
    st  = Vector{Float64}(undef, 2n - 1)
    slo = Vector{Float64}(undef, 2n - 1)
    shi = Vector{Float64}(undef, 2n - 1)
    for i in 1:n-1
        st[2i-1]  = t[i];   slo[2i-1] = lo[i]; shi[2i-1] = hi[i]
        st[2i]    = t[i+1]; slo[2i]   = lo[i]; shi[2i]   = hi[i]
    end
    st[2n-1] = t[n]; slo[2n-1] = lo[n]; shi[2n-1] = hi[n]
    st, slo, shi
end

# ── risktable_at ────────────────────────────────────────────────────────────

# Compute at-risk counts at specified time points from a KaplanMeier fit.
function risktable_at(km::KaplanMeier, timepoints)
    map(timepoints) do t
        idx = searchsortedfirst(km.events.time, t)
        idx <= length(km.events.natrisk) ? km.events.natrisk[idx] : 0
    end
end

# ── Kaplan-Meier transformation ────────────────────────────────────────────

Base.@kwdef struct KaplanMeierAnalysis
    interval::Union{Symbol,Nothing} = :confidence
    level::Float64 = 0.95
end

function (a::KaplanMeierAnalysis)(input::ProcessedLayer)
    all_times = reduce(vcat, input.positional[1])
    global_tmax = maximum(skipmissing(all_times))
    global_tp, _ = Makie.get_ticks(Makie.automatic, identity, Makie.automatic, 0.0, global_tmax)

    output = AlgebraOfGraphics.map(input) do p, n
        time, event = p
        km = fit(KaplanMeier, time, BitVector(event))
        ci = confint(km; level=a.level)

        t = [0.0; km.events.time]
        s = [1.0; km.survival]
        lo = [1.0; [c[1] for c in ci]]
        hi = [1.0; [c[2] for c in ci]]
        st, slo, shi = step_band_coords(t, lo, hi)
        nrisk = risktable_at(km, global_tp)

        return (t, s, st, slo, shi, global_tp, nrisk), (;)
    end

    ylabel = "Survival Probability"

    function set_ylabel(pl::ProcessedLayer)
        labels = copy(pl.labels)
        haskey(labels, 2) ? (labels[2] = ylabel) : insert!(labels, 2, ylabel)
        ProcessedLayer(pl; labels)
    end

    stairslayer = set_ylabel(ProcessedLayer(
        AlgebraOfGraphics.map(output) do p, n
            t, s, st, slo, shi, tp, nrisk = p
            (t, s), (;)
        end;
        plottype=Stairs, label=:survival,
        attributes=dictionary([:step => :post]),
    ))

    zerolayer = set_ylabel(ProcessedLayer(
        AlgebraOfGraphics.map(output) do p, n
            t, s, st, slo, shi, tp, nrisk = p
            ([first(t)], [0.0]), (;)
        end;
        plottype=Scatter, label=nothing,
        attributes=dictionary([:markersize => 0, :legend => (; visible=false)]),
    ))

    if isnothing(a.interval)
        return ProcessedLayers([stairslayer, zerolayer])
    end

    bandlayer = set_ylabel(ProcessedLayer(
        AlgebraOfGraphics.map(output) do p, n
            t, s, st, slo, shi, tp, nrisk = p
            (st, slo, shi), (;)
        end;
        plottype=Band, label=:ci,
        attributes=dictionary([:alpha => 0.3]),
    ))

    risktable_layer = set_ylabel(ProcessedLayer(
        AlgebraOfGraphics.map(output) do p, n
            t, s, st, slo, shi, tp, nrisk = p
            counts_str = [verbatim(join(string.(nrisk), ",")) for _ in tp]
            (tp, fill(NaN, length(tp))), (; inspector_label=counts_str)
        end;
        plottype=Scatter, label=nothing,
        attributes=dictionary([:visible => false, :legend => (; visible=false)]),
    ))

    return ProcessedLayers([bandlayer, stairslayer, risktable_layer, zerolayer])
end

"""
    kaplanmeier(; interval=:confidence, level=0.95)

Compute a Kaplan-Meier survival estimate as an AlgebraOfGraphics transformation.

The first positional mapping should be time, the second the event indicator
(0/1 or Bool).

Use `interval=:confidence` (default) for confidence bands, or
`interval=nothing` for just the survival curve. The `level` keyword
controls the confidence level (default 0.95).

Grouping (e.g. `color=:treatment`) is handled automatically by AlgebraOfGraphics.
"""
kaplanmeier(; options...) = transformation(KaplanMeierAnalysis(; options...))

# ── Censor ticks transformation ────────────────────────────────────────────

struct CensorTicksAnalysis end

function (::CensorTicksAnalysis)(input::ProcessedLayer)
    output = AlgebraOfGraphics.map(input) do p, n
        time, event = p
        km = fit(KaplanMeier, time, BitVector(event))
        cmask = km.events.ncensored .> 0
        ct = km.events.time[cmask]
        cs = km.survival[cmask]
        return (ct, cs), (;)
    end

    labels = copy(output.labels)
    haskey(labels, 2) ? (labels[2] = "Survival Probability") : insert!(labels, 2, "Survival Probability")

    return ProcessedLayer(output;
        plottype=Scatter,
        labels,
        attributes=dictionary([:markersize => 8]),
    )
end

"""
    censorticks()

Mark censored observations on a Kaplan-Meier plot. Use as a separate layer
combined with `kaplanmeier()`:

    plt = data(df) * mapping(:time, :event, color=:group, marker=:group) *
        (kaplanmeier() + censorticks())
"""
censorticks() = transformation(CensorTicksAnalysis())

# ── add_risktable! ──────────────────────────────────────────────────────────

function _extract_strata(ax)
    strata = Tuple{String,Vector{Float64},Vector{Int},Any}[]
    for plot in ax.scene.plots
        if plot isa Scatter && !to_value(plot.visible)
            pts = to_value(plot[1])
            tp = [p[1] for p in pts]
            label_val = to_value(plot.inspector_label)
            counts_str = label_val isa AbstractVector ? first(label_val) : label_val
            nrisk = parse.(Int, split(string(counts_str), ","))
            color = to_value(plot.color)
            push!(strata, ("", tp, nrisk, color))
        end
    end
    strata
end

function _legend_color_labels(fig)
    color_to_label = Dict{Any,String}()
    for c in fig.content
        c isa Makie.Legend || continue
        for (_, entries) in to_value(c.entrygroups)
            for entry in entries
                label = to_value(entry.label)
                for el in entry.elements
                    if el isa Makie.LineElement
                        color_to_label[to_value(el.linecolor)] = label
                    end
                end
            end
        end
    end
    color_to_label
end

"""
    add_risktable!(fg)

Add a risk table below a Kaplan-Meier plot. Call this on the result
returned by `draw()`. Works with both single and faceted layouts.
"""
function add_risktable!(fg)
    fig = fg.figure
    axes = filter(x -> x isa Axis, fig.content)
    isempty(axes) && return fg

    color_to_label = _legend_color_labels(fig)

    # Determine the grid row below the axes (assume axes are in row 1)
    ax_row = 2

    for (col_idx, ax) in enumerate(axes)
        strata = _extract_strata(ax)
        isempty(strata) && continue

        # Assign group names from legend
        for (i, (_, tp, nrisk, color)) in enumerate(strata)
            name = get(color_to_label, color, "")
            strata[i] = (name, tp, nrisk, color)
        end

        tp = strata[1][2]

        hidexdecorations!(ax, ticks=false, ticklabels=false, grid=false, minorgrid=false, label=false)

        nrows = length(strata) == 1 ? 2 : 1 + length(strata)
        show_labels = col_idx == 1
        row_labels = if length(strata) == 1
            show_labels ? ["-", "At Risk"] : [" ", " "]
        else
            if show_labels
                vcat([s[1] for s in reverse(strata)], ["At Risk"])
            else
                fill(" ", nrows)
            end
        end

        tab_ax = Axis(fig[ax_row, col_idx];
            xlabelvisible=false,
            height=nrows * 18,
            yticks=(1:nrows, row_labels),
            yticklabelsize=10, xticklabelsize=10,
            ylabelvisible=false,
            topspinevisible=false, bottomspinevisible=false,
            leftspinevisible=false, rightspinevisible=false,
            xgridvisible=false, ygridvisible=false,
            xticksvisible=false, xticklabelsvisible=false,
        )
        linkxaxes!(ax, tab_ax)

        for (j, (label, stp, nrisk, color)) in enumerate(reverse(strata))
            row = j
            c = length(strata) == 1 ? :black : color
            for (k, t) in enumerate(stp)
                text!(tab_ax, t, row; text=string(nrisk[k]),
                    align=(:center, :center), fontsize=9, color=c)
            end
        end

        ylims!(tab_ax, 0.5, nrows + 0.5)
    end

    fg
end

end # module
