"""
    plot_matrix_variable(var_data::Matrix, I::Vector{Int}, J::Vector{Int};
                         m=nothing, n=nothing, x=nothing, xlabel=nothing, var_name="",
                         vis_threshold=nothing, significance_fn=default_significance,
                         symmetric=false) -> Figure

Visualize a matrix variable (given in COO format) across multiple problem instances.

`var_data` has shape `(n_nonzeros, n_instances)`. `I[k]`/`J[k]` are the COO indices.

**Grid dimensions** (when not filtering):
- `m`/`n` fix the grid size; otherwise `maximum(I)`/`maximum(J)` is used.
- When `symmetric = true`, the grid is forced square; only one of `m`/`n` is needed.

**Thresholding**: when `n_nonzeros > vis_threshold`, only the top-vis_threshold most significant
rows *and* columns are shown in a compressed grid. Row/column significance = max entry score
in that row/column, where entry scores come from `significance_fn` (default: 1-norm).

**Symmetric mode** (`symmetric = true`): the COO is assumed to contain only one triangle.
Off-diagonal entries are always mirrored so the full symmetric subgrid is displayed.
When filtering, only the most significant columns are selected; the same indices serve as rows.
"""
function plot_matrix_variable(var_data::Matrix, I::Vector{Int}, J::Vector{Int};
                               m=nothing, n=nothing, x=nothing, xlabel=nothing, var_name="",
                               vis_threshold=nothing, significance_fn=default_significance,
                               symmetric=false)
    n_nonzeros, n_instances = size(var_data)
    @assert length(I) == n_nonzeros "Length of I ($(length(I))) must equal number of rows in var_data ($n_nonzeros)"
    @assert length(J) == n_nonzeros "Length of J ($(length(J))) must equal number of rows in var_data ($n_nonzeros)"

    # Resolve full matrix dimensions (used when not filtering)
    if symmetric
        if !isnothing(m) && !isnothing(n)
            @assert m == n "Symmetric matrix must be square (m=$m, n=$n)"
            effective_m_full = m; effective_n_full = n
        elseif !isnothing(m)
            effective_m_full = m; effective_n_full = m
        elseif !isnothing(n)
            effective_m_full = n; effective_n_full = n
        else
            d = max(maximum(I), maximum(J))
            effective_m_full = d; effective_n_full = d
        end
    else
        effective_m_full = isnothing(m) ? maximum(I) : m
        effective_n_full = isnothing(n) ? maximum(J) : n
    end

    @assert all(1 .<= I .<= effective_m_full) "All row indices must be in [1, m=$effective_m_full]"
    @assert all(1 .<= J .<= effective_n_full) "All column indices must be in [1, n=$effective_n_full]"

    x_vals = isnothing(x) ? collect(1:n_instances) : x
    x_label = isnothing(x) ? "Instance" : (isnothing(xlabel) ? "Unknown Parameter" : xlabel)

    entry_scores = [significance_fn(var_data[k, :]) for k in 1:n_nonzeros]
    I_plot, J_plot, data_row_idx, sel_rows, sel_cols, filtered =
        select_matrix_entries(entry_scores, I, J, vis_threshold, symmetric)

    data_plot = var_data[data_row_idx, :]
    n_plot = length(I_plot)

    # Compressed grid when filtering; full m×n grid otherwise
    if filtered
        row_map = Dict(r => idx for (idx, r) in enumerate(sel_rows))
        col_map = Dict(c => idx for (idx, c) in enumerate(sel_cols))
        n_grid_rows, n_grid_cols = length(sel_rows), length(sel_cols)
        grid_pos = (i, j) -> (row_map[i], col_map[j])
    else
        n_grid_rows, n_grid_cols = effective_m_full, effective_n_full
        grid_pos = (i, j) -> (i, j)
    end

    fig = Figure(size=(320 * n_grid_cols, 260 * n_grid_rows))
    gl = fig[1, 1] = GridLayout(n_grid_rows, n_grid_cols)

    for k in 1:n_plot
        entry_label = isempty(var_name) ? "[$(I_plot[k]),$(J_plot[k])]" : "$(var_name)[$(I_plot[k]),$(J_plot[k])]"
        gr, gc = grid_pos(I_plot[k], J_plot[k])
        ax = Axis(gl[gr, gc]; title=entry_label, xlabel=x_label, ylabel="Value")
        scatter!(ax, x_vals, data_plot[k, :])
    end

    return fig
end

"""
    plot_matrix_variable(var_data_dict::Dict{String, <:AbstractMatrix},
                         I::Vector{Int}, J::Vector{Int};
                         m=nothing, n=nothing, x=nothing, xlabel=nothing, var_name="",
                         vis_threshold=nothing, significance_fn=default_significance,
                         symmetric=false) -> Figure

Visualize a matrix variable from multiple solvers in a COO grid layout.

Each key in `var_data_dict` is a solver name; its value is a `(n_nonzeros, n_instances)` matrix.
`I` and `J` are the shared sparsity structure. All other options behave as in the single-solver
method. Entry significance is computed from values concatenated across all solvers.
"""
function plot_matrix_variable(var_data_dict::Dict{String, <:AbstractMatrix},
                               I::Vector{Int}, J::Vector{Int};
                               m=nothing, n=nothing, x=nothing, xlabel=nothing, var_name="",
                               vis_threshold=nothing, significance_fn=default_significance,
                               symmetric=false)
    solver_names = collect(keys(var_data_dict))
    first_data = var_data_dict[first(solver_names)]
    n_nonzeros, n_instances = size(first_data)

    if symmetric
        if !isnothing(m) && !isnothing(n)
            @assert m == n "Symmetric matrix must be square (m=$m, n=$n)"
            effective_m_full = m; effective_n_full = n
        elseif !isnothing(m)
            effective_m_full = m; effective_n_full = m
        elseif !isnothing(n)
            effective_m_full = n; effective_n_full = n
        else
            d = max(maximum(I), maximum(J))
            effective_m_full = d; effective_n_full = d
        end
    else
        effective_m_full = isnothing(m) ? maximum(I) : m
        effective_n_full = isnothing(n) ? maximum(J) : n
    end

    @assert all(1 .<= I .<= effective_m_full) "All row indices must be in [1, m=$effective_m_full]"
    @assert all(1 .<= J .<= effective_n_full) "All column indices must be in [1, n=$effective_n_full]"

    x_vals = isnothing(x) ? collect(1:n_instances) : x
    x_label = isnothing(x) ? "Instance" : (isnothing(xlabel) ? "Unknown Parameter" : xlabel)

    # Combine all solver data for significance scoring
    entry_scores = [significance_fn(vcat([var_data_dict[s][k, :] for s in solver_names]...))
                    for k in 1:n_nonzeros]
    I_plot, J_plot, data_row_idx, sel_rows, sel_cols, filtered =
        select_matrix_entries(entry_scores, I, J, vis_threshold, symmetric)

    n_plot = length(I_plot)

    if filtered
        row_map = Dict(r => idx for (idx, r) in enumerate(sel_rows))
        col_map = Dict(c => idx for (idx, c) in enumerate(sel_cols))
        n_grid_rows, n_grid_cols = length(sel_rows), length(sel_cols)
        grid_pos = (i, j) -> (row_map[i], col_map[j])
    else
        n_grid_rows, n_grid_cols = effective_m_full, effective_n_full
        grid_pos = (i, j) -> (i, j)
    end

    palette = Makie.wong_colors()
    solver_colors = Dict(name => palette[mod1(i, length(palette))]
                         for (i, name) in enumerate(solver_names))

    fig = Figure(size=(320 * n_grid_cols, 260 * n_grid_rows + 60))
    gl = fig[1, 1] = GridLayout(n_grid_rows, n_grid_cols)

    legend_handles = []
    legend_labels = String[]

    for k in 1:n_plot
        entry_label = isempty(var_name) ? "[$(I_plot[k]),$(J_plot[k])]" : "$(var_name)[$(I_plot[k]),$(J_plot[k])]"
        gr, gc = grid_pos(I_plot[k], J_plot[k])
        ax = Axis(gl[gr, gc]; title=entry_label, xlabel=x_label, ylabel="Value")

        for solver_name in solver_names
            p = scatter!(ax, x_vals, var_data_dict[solver_name][data_row_idx[k], :];
                         color=solver_colors[solver_name])
            if k == 1
                push!(legend_handles, p)
                push!(legend_labels, solver_name)
            end
        end
    end

    Legend(fig[2, 1], legend_handles, legend_labels;
           orientation=:horizontal, tellwidth=false)

    return fig
end
