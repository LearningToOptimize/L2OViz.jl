"""
    plot_matrix_variable(I::Vector{Int}, J::Vector{Int}, x, var_data::Matrix...;
                         solver_names=nothing,
                         xlabel=nothing, var_name="",
                         vis_threshold::Int=20, significance_fn=default_significance) -> Figure

Visualize a symmetric matrix variable (given in COO format) across multiple problem instances.
The matrix variable is assumed to have the same COO across all the problem instances.
The COO coordinates specified by `I` and `J` should not contain both (i, j) and (j, i), and entries are
visualized for only the given half of the matrix specified by the coordinates.

`x` is either:
- A single `Vector`: the same x is used for every solver's data matrix. Length must equal the
number of instances in the data of each solver.
- Multiple `Vector`s: one `Vector` per solver, in the same order as `var_data`. `length(x)`
must equal the number of solvers, and length of each Vector must equal the number of instances in
the data of each solver.

`var_data` is one or more `(nnz × n_instances)` data matrices, one per solver.

If not provided, solver names default to "Solver 1", "Solver 2", ...

The x-axis label defaults to `"Unknown Parameter"` unless `xlabel` is given.

**Thresholding**: if number of unique rows/columns exceeds `vis_threshold`, select an induced
submatrix with dimension `vis_threshold`. `significance_fn` (default: 1-norm) is applied
to the values of the k-th entry of the variable across all solvers and instances to get the
score of coordinate `(I[k], J[k])`, and the score of each column/row is the max entry score in
that column/row. If there are repeated `(i, j)` coordinates, only the highest-scoring one is
kept after thresholding (with a warning).

**Grid dimensions**: the side length of the (square) matrix is taken as
`max(maximum(I), maximum(J))`.
"""
function plot_matrix_variable(I::Vector{Int}, J::Vector{Int}, x, var_data::Matrix...;
                              solver_names=nothing,
                              xlabel=nothing, var_name="",
                              vis_threshold::Int=20,
                              significance_fn=default_significance)
    length(var_data) >= 1 || throw(ArgumentError("At least one data matrix must be provided"))
    nnz = validate_var_data_dims(var_data)
    length(I) == nnz || throw(DimensionMismatch("Length of I must equal number of rows in var_data"))
    length(J) == nnz || throw(DimensionMismatch("Length of J must equal number of rows in var_data"))
    x_vecs = resolve_x_vecs(x, var_data)
    n_solvers = length(var_data)
    solver_names = resolve_solver_names(solver_names, n_solvers)
    x_label = isnothing(xlabel) ? "Unknown Parameter" : xlabel

    # Full matrix dimension; symmetric matrix is square so the side length is the
    # largest index appearing in either coordinate.
    n = max(maximum(I), maximum(J))

    # For each entry, combine values across all solvers and all instances for significance score
    entry_scores = compute_entry_scores(var_data, nnz, significance_fn)
    I_plot, J_plot, selected_indices, filtered =
        select_matrix_entries(entry_scores, I, J, vis_threshold)

    n, grid_pos = matrix_grid_layout(I_plot, J_plot, filtered, n)

    solver_colors = solver_palette(n_solvers)

    fig = Figure(size=(320 * n, 260 * n + 60))
    gl = fig[1, 1] = GridLayout(n, n)

    axes, legend_handles = draw_matrix_panels!(
        gl, var_data, x_vecs, x_label, var_name,
        I_plot, J_plot, selected_indices, grid_pos, solver_colors)

    linkxaxes!(axes...)
    linkyaxes!(axes...)

    Legend(fig[2, 1], legend_handles, solver_names;
           orientation=:horizontal, tellwidth=false)

    return fig
end
