"""
    plot_polar_variable(x, var_data::Tuple{AbstractMatrix, AbstractMatrix}...;
                        solver_names=nothing, xlabel=nothing, var_name="",
                        vis_threshold::Int=20, significance_fn=default_significance,
                        symlog::Bool=false, palette=nothing, alpha=1.0) -> Figure

Visualize a complex-valued vector variable across multiple problem instances in polar
coordinates. Each entry is plotted on its own polar axis, with radius `|z|` and angle `arg(z)`.

`x` is either:
- A single `Vector`: the same x is used for every solver's data. Length must equal the number of
instances in the data of each solver.
- Multiple `Vector`s: one `Vector` per solver, in the same order as `var_data`. `length(x)`
must equal the number of solvers, and length of each Vector must equal the number of instances in
the data of each solver.

`x` is encoded as the marker size, which grows linearly with `x` over the range of `x` across all
solvers.

`var_data` is one or more `(real_part, imag_part)` tuples, one per solver, where both parts are
`(n_entries × n_instances)` matrices.

If not provided, solver names default to "Solver 1", "Solver 2", ...

The marker-size legend title defaults to `"Unknown Parameter"` unless `xlabel` is given.

Only the top-`vis_threshold` most significant entries are visualized. `significance_fn`
(default: 1-norm) is applied to the radii of each entry across all solvers and instances to get
the score of each entry. Entry labels always reflect the original indices.

Set `symlog=true` to draw the radii on a symmetric log scale, which is logarithmic above the
smallest nonzero displayed radius and keeps zero radii at the origin.

`palette` sets the per-solver colors; see [`plot_variable`](@ref).
"""
function plot_polar_variable(x, var_data::Tuple{AbstractMatrix, AbstractMatrix}...;
                             solver_names=nothing, xlabel=nothing,
                             var_name="", vis_threshold::Int=20,
                             significance_fn=default_significance,
                             symlog::Bool=false, palette=nothing, alpha::Real=1.0)
    radius_data, angle_data, n_entries, x_vecs, solver_names, x_label =
        prepare_polar_inputs(x, var_data, solver_names, xlabel)

    # For each entry, combine radii across all solvers and all instances for significance score
    entry_scores = compute_entry_scores(radius_data, n_entries, significance_fn)
    selected_indices = select_variable_entries(entry_scores, vis_threshold)

    n_rows, n_cols = vector_grid_layout(length(selected_indices))
    panel_cells = [vector_grid_cell(k, n_cols) for k in eachindex(selected_indices)]
    panel_titles = [entry_label(var_name, data_idx) for data_idx in selected_indices]

    return polar_figure(n_rows, n_cols, radius_data, angle_data, x_vecs, panel_cells,
                        panel_titles, selected_indices, solver_names, x_label, palette;
                        symlog=symlog, alpha=alpha)
end

"""
    plot_polar_graph_variable(I::Vector{Int}, J::Vector{Int}, x,
                              var_data::Tuple{AbstractMatrix, AbstractMatrix}...;
                              solver_names=nothing, xlabel=nothing, var_name="",
                              vis_threshold::Int=20, significance_fn=default_significance,
                              symlog::Bool=false, palette=nothing, alpha=1.0, bipartite=false,
                              scored_on=:rows) -> Figure

Visualize a complex-valued graph variable across multiple problem instances in polar
coordinates. Edge `k`, given by the vertex pair `(I[k], J[k])`, is plotted on a polar axis in grid
cell `(I[k], J[k])`, with radius `|z|` and angle `arg(z)`. The grid layout, `bipartite` and
`scored_on` follow [`plot_graph_variable`](@ref).

`x`, `var_data`, `solver_names`, `xlabel`, `symlog` and `palette` are as in
[`plot_polar_variable`](@ref), with `n_e` (the number of edges) rows in each data matrix.

**Thresholding**: as in [`plot_graph_variable`](@ref), with `significance_fn` (default: 1-norm)
applied to the radii of each entry across all solvers and instances.
"""
function plot_polar_graph_variable(I::Vector{Int}, J::Vector{Int}, x,
                                   var_data::Tuple{AbstractMatrix, AbstractMatrix}...;
                                   solver_names=nothing, xlabel=nothing,
                                   var_name="", vis_threshold::Int=20,
                                   significance_fn=default_significance,
                                   symlog::Bool=false, palette=nothing, alpha::Real=1.0,
                                   bipartite::Bool=false, scored_on::Symbol=:rows)
    radius_data, angle_data, n_e, x_vecs, solver_names, x_label =
        prepare_polar_inputs(x, var_data, solver_names, xlabel)
    length(I) == n_e || throw(DimensionMismatch("Length of I must equal number of rows in var_data"))
    length(J) == n_e || throw(DimensionMismatch("Length of J must equal number of rows in var_data"))

    # For each entry, combine radii across all solvers and all instances for significance score
    entry_scores = compute_entry_scores(radius_data, n_e, significance_fn)
    I_plot, J_plot, selected_indices = select_variable_edges(entry_scores, I, J, vis_threshold;
                                                               bipartite=bipartite,
                                                               scored_on=scored_on)

    n_rows, n_cols, grid_pos = matrix_grid_layout(I_plot, J_plot; bipartite=bipartite)
    panel_cells = [grid_pos(i, j) for (i, j) in zip(I_plot, J_plot)]
    panel_titles = [entry_label(var_name, i, j) for (i, j) in zip(I_plot, J_plot)]

    return polar_figure(n_rows, n_cols, radius_data, angle_data, x_vecs, panel_cells,
                        panel_titles, selected_indices, solver_names, x_label, palette;
                        symlog=symlog, alpha=alpha)
end

# Validate the (real, imag) tuples and convert them to polar coordinates. Returns per-solver
# radius and angle matrices, the variable dimension, resolved x vectors, solver names and x label.
function prepare_polar_inputs(x, var_data, solver_names, xlabel)
    length(var_data) >= 1 || throw(ArgumentError("At least one data tuple must be provided"))
    validate_polar_pair_dims(var_data)
    polar_data = [polar_coordinates(real_part, imag_part) for (real_part, imag_part) in var_data]
    radius_data = first.(polar_data)
    angle_data = last.(polar_data)
    n_entries = validate_var_data_dims(radius_data)
    x_vecs = resolve_x_vecs(x, radius_data)
    solver_names = resolve_solver_names(solver_names, length(var_data))
    x_label = isnothing(xlabel) ? "Unknown Parameter" : xlabel
    return radius_data, angle_data, n_entries, x_vecs, solver_names, x_label
end

# Build a figure with an `n_rows × n_cols` grid of polar panels and the legend below it.
function polar_figure(n_rows, n_cols, radius_data, angle_data, x_vecs, panel_cells,
                      panel_titles, data_indices, solver_names, x_label, palette;
                      symlog::Bool=false, alpha=1.0)
    solver_colors = solver_palette(length(radius_data), palette)
    marker_sizes, size_key = x_marker_sizes(x_vecs)
    rscale = resolve_yscale(symlog, radius_data, data_indices)

    # Shared radius range over the displayed entries with 5% headroom (in transformed space) so
    # markers at the largest radius are not clipped by the rim; fall back to 1 if every radius
    # is zero.
    _, max_radius = data_extrema(radius_data, data_indices)
    rmax_scaled = max_radius > 0 ? 1.05 * rscale(max_radius) : 1.0

    fig = Figure(size=(320 * n_cols, 320 * n_rows + 100))
    panel_grid = fig[1, 1] = GridLayout(n_rows, n_cols)
    draw_polar_panels!(panel_grid, radius_data, angle_data, marker_sizes, panel_cells,
                       panel_titles, data_indices, solver_colors, rmax_scaled;
                       rscale=rscale, alpha=alpha)
    add_polar_legend!(fig[2, 1], solver_colors, solver_names, size_key, x_label)
    return fig
end
