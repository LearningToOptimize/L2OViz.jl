# Coordinates and values of the nonzero entries of `d`, with `x_vec[j]` as the x-coordinate of
# instance `j` and the entry index as the y-coordinate. Entries with `abs(value) <= atol` count
# as zero.
function nonzero_points(d::AbstractMatrix, x_vec::AbstractVector, atol::Real)
    nonzero_coords = findall(value -> abs(value) > atol, d)
    xs = [x_vec[c[2]] for c in nonzero_coords]
    entry_indices = [c[1] for c in nonzero_coords]
    values = [d[c] for c in nonzero_coords]
    return xs, entry_indices, values
end

# Shared color range over every solver's displayed values. A degenerate range (a single
# distinct value, or no values at all) is padded so the color bar stays well defined.
function sparse_color_range(all_values)
    isempty(all_values) && return (0.0, 1.0)
    lo, hi = extrema(all_values)
    if lo == hi
        pad = lo == 0 ? 0.5 : 0.5 * abs(lo)
        return (lo - pad, hi + pad)
    end
    return (lo, hi)
end

# Integer y-ticks for entry indices 1:n_entries, at a 1/2/5 × 10^k step giving at most
# `max_ticks` ticks.
function entry_index_ticks(n_entries::Int; max_ticks::Int=10)
    magnitude = 10^max(0, floor(Int, log10(n_entries / max_ticks)))
    step = first(m * magnitude for m in (1, 2, 5, 10) if n_entries / (m * magnitude) <= max_ticks)
    first_tick = step == 1 ? 1 : step
    return first_tick:step:n_entries
end

# Color bar ticks. Makie's automatic ticks leave the negative half of a symlog bar unlabeled,
# so for symlog we place ticks at zero and at ±10^k for every decade from the linear threshold
# up to the largest magnitude within `colorrange`.
sparse_colorbar_ticks(::Any, colorrange) = Makie.automatic
function sparse_colorbar_ticks(colorscale::Makie.Symlog10, colorrange)
    lo, hi = colorrange
    max_abs = max(abs(lo), abs(hi))
    # Symlog10 stores the linear threshold as `lower == -upper`.
    first_decade = ceil(Int, log10(colorscale.upper))
    last_decade = floor(Int, log10(max_abs))
    decades = first_decade:last_decade
    signed_powers = [(sign, k) for sign in (-1, 1) for k in decades if lo <= sign * 10.0^k <= hi]
    sort!(signed_powers; by=((sign, k),) -> sign * 10.0^k)
    tick_values = [sign * 10.0^k for (sign, k) in signed_powers]
    tick_labels = Any[rich(sign < 0 ? "−10" : "10", superscript(string(k)))
                      for (sign, k) in signed_powers]
    if lo <= 0 <= hi
        # Insert 0 between the negative and positive ticks.
        zero_position = searchsortedfirst(tick_values, 0.0)
        insert!(tick_values, zero_position, 0.0)
        insert!(tick_labels, zero_position, "0")
    end
    return (tick_values, tick_labels)
end

"""
    plot_sparse_variable(x, var_data::AbstractMatrix...; solver_names=nothing, xlabel=nothing,
                         var_name="", symlog::Bool=false, atol::Real=1e-6, colormap=:viridis,
                         markersize::Real=6, alpha=1.0) -> Figure

Visualize the sparsity pattern and values of a vector variable across multiple problem
instances. Every nonzero entry is drawn as a scatter point at (instance `x`, entry index),
colored by its value; zero entries are omitted.

`x` is either:
- A single `Vector`: the same x is used for every solver's data matrix. Length must equal the
number of instances in the data of each solver.
- Multiple `Vector`s: one `Vector` per solver, in the same order as `var_data`. `length(x)`
must equal the number of solvers, and length of each Vector must equal the number of instances in
the data of each solver.

`var_data` is one or more `(n_entries × n_instances)` matrices, one per solver. Each solver gets
its own panel; the panels are placed side-by-side and share the y-axis (entry index) and a single
color bar.

If not provided, solver names default to "Solver 1", "Solver 2", ...

The x-axis label defaults to `"Unknown Parameter"` unless `xlabel` is given.

Only entries with `abs(value) > atol` are displayed (default `atol=1e-6`).

Set `symlog=true` to draw the color bar on a symmetric log scale, which behaves as a log scale
for magnitudes above the smallest displayed magnitude while still supporting negative values.

`colormap` is any colormap Makie accepts (e.g. `:viridis`, `:balance`).
"""
function plot_sparse_variable(x, var_data::AbstractMatrix...;
                              solver_names=nothing, xlabel=nothing, var_name="",
                              symlog::Bool=false, atol::Real=1e-6, colormap=:viridis,
                              markersize::Real=6, alpha::Real=1.0)
    length(var_data) >= 1 || throw(ArgumentError("At least one data matrix must be provided"))
    atol >= 0 || throw(ArgumentError("atol must be nonnegative"))
    n_entries = validate_var_data_dims(var_data)
    x_vecs = resolve_x_vecs(x, var_data)
    n_solvers = length(var_data)
    solver_names = resolve_solver_names(solver_names, n_solvers)
    x_label = isnothing(xlabel) ? "Unknown Parameter" : xlabel

    solver_points = [nonzero_points(d, x_vecs[i], atol) for (i, d) in enumerate(var_data)]
    all_values = reduce(vcat, (points[3] for points in solver_points))
    isempty(all_values) && @warn "plot_sparse_variable: every entry is zero; panels will be empty."

    colorscale = resolve_symlog_transform(symlog && !isempty(all_values), all_values)
    colorrange = sparse_color_range(all_values)

    # Panels are taller than the other scatter plots, since each one condenses every entry.
    panel_width = 450
    panel_height = 650
    fig = Figure(size=(panel_width * n_solvers + 150, panel_height + 80))

    y_label = isempty(var_name) ? "Entry Index" : "$(var_name) Entry Index"
    axes = Axis[]
    scatter_plot = nothing
    for (i, (xs, entry_indices, values)) in enumerate(solver_points)
        ax = Axis(fig[1, i]; title=solver_names[i], xlabel=x_label,
                  ylabel=(i == 1 ? y_label : ""), yticks=entry_index_ticks(n_entries))
        i > 1 && hideydecorations!(ax; grid=false, ticks=false)
        push!(axes, ax)
        scatter_plot = scatter!(ax, xs, entry_indices; color=values, colormap=colormap,
                                colorrange=colorrange, colorscale=colorscale,
                                markersize=markersize, alpha=alpha)
    end

    linkxaxes!(axes...)
    linkyaxes!(axes...)
    ylims!(axes[1], 0.5, n_entries + 0.5)

    color_label = isempty(var_name) ? "Value" : var_name
    Colorbar(fig[1, n_solvers + 1], scatter_plot; label=color_label,
             ticks=sparse_colorbar_ticks(colorscale, colorrange))
    return fig
end
