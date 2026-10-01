# Normalize a single positional variable argument into a Vector of its per-solver matrices.
# A lone Matrix is a single-solver variable; a Vector/Tuple of Matrices lists a variable's
# solvers in shared solver order.
variable_solver_matrices(d::AbstractMatrix) = [d]
function variable_solver_matrices(group)
    if group isa Union{AbstractVector, Tuple} && all(g -> g isa AbstractMatrix, group)
        return collect(group)
    end
    throw(ArgumentError(
        "Each variable must be a Matrix (single solver) or a Vector/Tuple of Matrices " *
        "(multiple solvers); got $(typeof(group))"))
end

# Per-entry geometric mean of the absolute values across instances (the second dimension),
# computed in log space. An entry with any zero has geometric mean zero (`log(0) == -Inf`
# propagates to `exp(-Inf) == 0`), which is the correct limit.
geomean_abs(d::AbstractMatrix) = exp.(vec(mean(log.(abs.(d)), dims=2)))

"""
    plot_scale_histograms(var_data...; solver_names=nothing, var_names=nothing,
                          xlabel=nothing, symlog=false, bins::Int=20,
                          palette=nothing, alpha=1.0) -> Figure

Visualize the distribution of per-entry scales of one or more variables, comparing solvers.

Each positional argument is one variable, given as either:
- A single `Matrix`: the variable has a single solver.
- A `Vector`/`Tuple` of `Matrix`: the variable's per-solver data, in the same solver order as
`solver_names`.

Every solver matrix has shape `(y_dim × n_instances)`. The matrices of the same variable must
share dimensions across solvers, but different variables may differ in dimensions. All
variables must list the same number of solvers.

For each solver matrix, every entry (row) is summarized across instances by four statistics,
each a length-`y_dim` vector: the mean of the absolute values, the geometric mean of the
absolute values, the mean of the original values, and the median of the original values. The
output is a grid with one row per variable and four columns in that order; each panel is a
histogram of the corresponding per-entry values. Within a panel the solvers' histograms are
drawn as side-by-side (dodged) bars. The two absolute-value panels of a variable share a common
x-axis, and the two original-value panels share another.

If not provided, solver names default to "Solver 1", "Solver 2", ..., and variable names to
"Variable 1", "Variable 2", ....

The x-axis label defaults to `"Variable Scale"` unless `xlabel` is given.

Set `symlog=true` to draw the x-axis (variable scale) on a symmetric log scale for every
variable; the histogram bins are then spaced evenly in the transformed space so the bars keep an
even visual width. Pass a `Vector{Bool}` of length equal to the number of variables to choose
the scale per variable.

`palette` sets the per-solver colors; it defaults to `Makie.wong_colors()`. It may be a vector
of colors (of any type Makie accepts, e.g. `[:red, :blue]`), which is cycled through when there
are more solvers than colors, or a `Symbol` naming a Makie/ColorSchemes palette (e.g. `:tab10`,
`:viridis`): categorical palettes use their discrete colors, and continuous colormaps are sampled
into as many evenly spaced colors as there are solvers.
"""
function plot_scale_histograms(var_data...; solver_names=nothing, var_names=nothing,
                               xlabel=nothing, symlog=false, bins::Int=20,
                               palette=nothing, alpha::Real=1.0)
    length(var_data) >= 1 || throw(ArgumentError("At least one variable must be provided"))
    bins > 0 || throw(ArgumentError("bins must be positive"))

    # One entry per variable, each a Vector of that variable's per-solver matrices.
    variables = [variable_solver_matrices(group) for group in var_data]
    n_vars = length(variables)
    n_solvers = length(variables[1])
    for (v, solvers) in enumerate(variables)
        length(solvers) == n_solvers || throw(ArgumentError(
            "Variable $v lists $(length(solvers)) solver(s); expected $n_solvers to match variable 1"))
        # Matrices of the same variable must share dimensions across solvers.
        validate_var_data_dims(solvers)
        for (i, d) in enumerate(solvers)
            size(d, 2) >= 1 || throw(ArgumentError(
                "Variable $v, solver $i has no instances (size(d, 2) == 0)"))
        end
    end

    solver_names = resolve_names(solver_names, n_solvers, "solver_names", "Solver")
    var_names = resolve_names(var_names, n_vars, "var_names", "Variable")
    x_label = isnothing(xlabel) ? "Variable Scale" : xlabel
    solver_colors = solver_palette(n_solvers, palette)
    symlog_flags = resolve_symlog_flags(symlog, n_vars)

    label_col = 1
    stat_headers = ("Mean of abs", "Geomean", "Mean", "Median")
    stat_fns = (d -> vec(mean(abs.(d), dims=2)),
                geomean_abs,
                d -> vec(mean(d, dims=2)),
                d -> vec(median(d, dims=2)))
    # The first two absolute-value statistics share one scale; the remaining two share another.
    stat_groups = ((1, 2), (3, 4))
    stat_cols = (1:length(stat_headers)) .+ label_col
    first_stat_col, last_stat_col = first(stat_cols), last(stat_cols)

    fig = Figure(size=(360 * length(stat_cols) + 80, 260 * n_vars + 120))

    for (col, header) in zip(stat_cols, stat_headers)
        Label(fig[1, col], header; font=:bold, tellwidth=false)
    end

    for (v, solvers) in enumerate(variables)
        grid_row = v + 1
        Label(fig[grid_row, label_col], var_names[v];
              font=:bold, rotation=pi / 2, tellheight=false)
        # Per-statistic, per-solver per-entry summaries (each with length y_dim)
        stat_series = [[stat_fn(d) for d in solvers] for stat_fn in stat_fns]
        for group in stat_groups
            group_series = [stat_series[i] for i in group]
            all_values = reduce(vcat, (reduce(vcat, series) for series in group_series))
            transform = resolve_symlog_transform(symlog_flags[v], all_values)
            edges = scale_histogram_bins(all_values, bins, transform)
            group_axes = Axis[]
            for i in group
                col = stat_cols[i]
                ax = Axis(fig[grid_row, col]; xscale=transform, xlabel=x_label,
                          ylabel=(col == first_stat_col ? "Count" : ""))
                hist!(ax, stat_series[i]; bins=edges, dodge=1:n_solvers, n_dodge=n_solvers,
                      color=solver_colors, alpha=alpha)
                push!(group_axes, ax)
            end
            linkxaxes!(group_axes...)
        end
    end

    legend_handles = [PolyElement(; color=solver_colors[i]) for i in 1:n_solvers]
    Legend(fig[n_vars + 2, first_stat_col:last_stat_col], legend_handles, solver_names;
           orientation=:horizontal, tellwidth=false)

    return fig
end
