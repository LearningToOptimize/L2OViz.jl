"""
    plot_variable(var_data::Matrix; x=nothing, xlabel=nothing, var_name="",
                  vis_threshold=nothing, significance_fn=default_significance) -> Figure

Visualize each entry of a vector variable across multiple problem instances.

`var_data` has shape `(n_entries, n_instances)`. Each row becomes one scatter subplot tiled
into a roughly square grid. The x-axis defaults to instance indices with label `"Instance"`;
if custom `x` values are provided the label defaults to `"Unknown Parameter"` unless `xlabel`
is also given.

If `vis_threshold` is provided and `n_entries > vis_threshold`, only the `vis_threshold` entries
with the highest significance scores (via `significance_fn`) are shown. Entry labels always
reflect original indices.
"""
function plot_variable(var_data::Matrix; x=nothing, xlabel=nothing, var_name="",
                       vis_threshold=nothing, significance_fn=default_significance)
    n_entries, n_instances = size(var_data)
    x_vals = isnothing(x) ? collect(1:n_instances) : x
    x_label = isnothing(x) ? "Instance" : (isnothing(xlabel) ? "Unknown Parameter" : xlabel)

    if !isnothing(vis_threshold) && n_entries > vis_threshold
        scores = [significance_fn(var_data[k, :]) for k in 1:n_entries]
        selected_indices, _ = select_variable_entries(scores, vis_threshold)
    else
        selected_indices = collect(1:n_entries)
    end

    n_plot = length(selected_indices)
    n_cols = ceil(Int, sqrt(n_plot))
    n_rows = ceil(Int, n_plot / n_cols)

    fig = Figure(size=(320 * n_cols, 260 * n_rows))

    for (plot_k, orig_k) in enumerate(selected_indices)
        grid_row = div(plot_k - 1, n_cols) + 1
        grid_col = mod(plot_k - 1, n_cols) + 1
        entry_label = isempty(var_name) ? "[$orig_k]" : "$(var_name)[$orig_k]"
        ax = Axis(fig[grid_row, grid_col]; title=entry_label, xlabel=x_label, ylabel="Value")
        scatter!(ax, x_vals, var_data[orig_k, :])
    end

    return fig
end

"""
    plot_variable(var_data_dict::Dict{String, <:AbstractMatrix}; x=nothing, xlabel=nothing,
                  var_name="", vis_threshold=nothing, significance_fn=default_significance) -> Figure

Visualize a vector variable from multiple solvers on the same set of subplots.

Each key in `var_data_dict` is a solver name; its value is a `(n_entries, n_instances)` matrix.
All solver data is overlaid on each subplot using distinct colors. A shared legend is appended
below the subplot grid. The x-axis label follows the same rules as the single-solver method.

Entry significance is computed by applying `significance_fn` to each entry's values concatenated
across all solvers.
"""
function plot_variable(var_data_dict::Dict{String, <:AbstractMatrix}; x=nothing, xlabel=nothing,
                       var_name="", vis_threshold=nothing, significance_fn=default_significance)
    solver_names = collect(keys(var_data_dict))
    first_data = var_data_dict[first(solver_names)]
    n_entries, n_instances = size(first_data)
    x_vals = isnothing(x) ? collect(1:n_instances) : x
    x_label = isnothing(x) ? "Instance" : (isnothing(xlabel) ? "Unknown Parameter" : xlabel)

    if !isnothing(vis_threshold) && n_entries > vis_threshold
        scores = [significance_fn(vcat([var_data_dict[s][k, :] for s in solver_names]...))
                  for k in 1:n_entries]
        selected_indices, _ = select_variable_entries(scores, vis_threshold)
    else
        selected_indices = collect(1:n_entries)
    end

    n_plot = length(selected_indices)
    n_cols = ceil(Int, sqrt(n_plot))
    n_rows = ceil(Int, n_plot / n_cols)

    palette = Makie.wong_colors()
    # Map each solver to a fixed color
    solver_colors = Dict(name => palette[mod1(i, length(palette))]
                         for (i, name) in enumerate(solver_names))

    # Extra height for the legend row
    fig = Figure(size=(320 * n_cols, 260 * n_rows + 60))

    legend_handles = []
    legend_labels = String[]

    for (plot_k, orig_k) in enumerate(selected_indices)
        grid_row = div(plot_k - 1, n_cols) + 1
        grid_col = mod(plot_k - 1, n_cols) + 1
        entry_label = isempty(var_name) ? "[$orig_k]" : "$(var_name)[$orig_k]"
        ax = Axis(fig[grid_row, grid_col]; title=entry_label, xlabel=x_label, ylabel="Value")

        for solver_name in solver_names
            p = scatter!(ax, x_vals, var_data_dict[solver_name][orig_k, :];
                         color=solver_colors[solver_name])
            # Collect handles from the first subplot for the shared legend
            if plot_k == 1
                push!(legend_handles, p)
                push!(legend_labels, solver_name)
            end
        end
    end

    Legend(fig[n_rows + 1, 1:n_cols], legend_handles, legend_labels;
           orientation=:horizontal, tellwidth=false)

    return fig
end
