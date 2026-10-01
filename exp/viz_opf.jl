# exp/viz_opf.jl
# Utility function for visualizing OPF solution variables.

using CairoMakie
using L2OViz
using PGLib
using PowerModels

# Load branch and bus-pair connectivity (from-bus/to-bus indices) plus component counts
# from a PowerModels network dict. Bus pairs are the unique branch endpoint pairs, so
# parallel branches sharing endpoints collapse to a single pair.
function _get_power_system_data(network::Dict)
    branch_dict = network["branch"]
    sorted_branch_keys = sort(collect(keys(branch_dict)), by=k -> parse(Int, string(k)))
    I_branches = [branch_dict[k]["f_bus"] for k in sorted_branch_keys]
    J_branches = [branch_dict[k]["t_bus"] for k in sorted_branch_keys]
    n_branches = length(sorted_branch_keys)
    n_buses    = length(network["bus"])

    # Bus pairs are the unique (f_bus, t_bus) branch endpoints, kept in order of first
    # occurrence among the sorted branches; parallel branches collapse to a single pair.
    unique_buspairs = unique(zip(I_branches, J_branches))
    I_buspairs = [buspair[1] for buspair in unique_buspairs]
    J_buspairs = [buspair[2] for buspair in unique_buspairs]
    n_buspairs = length(unique_buspairs)

    return (; I_branches, J_branches, n_branches, n_buses, I_buspairs, J_buspairs, n_buspairs)
end

_load_network(system_identifier::String, flat::Bool) =
    flat ? Dict() : make_basic_network(pglib(system_identifier))

# Group the per-solver `var_data` by variable into `(var_name, per-solver arrays)` pairs.
# Array mode: `variables` is a String naming the single variable whose per-solver arrays are
# `var_data`.
# Dict mode: `variables` is a Vector{String} of keys into every per-solver Dict.
# Every array must have an ndims value in `allowed_ndims`.
function _var_data_pairs(variables, var_data, allowed_ndims)
    length(var_data) >= 1 || throw(ArgumentError("At least one var_data element is required"))
    all_dict  = all(d -> d isa Dict, var_data)
    all_array = all(d -> d isa AbstractArray, var_data)
    all_dict || all_array || throw(ArgumentError(
        "All var_data elements must be Dicts, or all must be arrays"))

    if all_array
        variables isa String || throw(ArgumentError(
            "`variables` should be a String when var_data are arrays (single variable)"))
        var_data_pairs = [(variables, collect(var_data))]
    else
        variables isa Vector{String} || throw(ArgumentError(
            "`variables` should be a Vector{String} when var_data are Dicts (multiple variables)"))
        # Check that all solvers have all variables.
        for (i, d) in enumerate(var_data)
            for v in variables
                haskey(d, v) || throw(ArgumentError("Variable '$v' not found in data for solver $(i)"))
            end
        end
        var_data_pairs = [(v, [d[v] for d in var_data]) for v in variables]
    end

    for (var_name, solvers_data) in var_data_pairs, (i, d) in enumerate(solvers_data)
        d isa AbstractArray && ndims(d) in allowed_ndims || throw(ArgumentError(
            "Data of variable '$var_name' for solver $i must be an array with ndims in " *
            "$allowed_ndims, got $(typeof(d))"))
    end
    return var_data_pairs
end

# Graph endpoints `(I, J)` of a variable with dimension `n_dim`:
# - branch endpoints for a branch variable
# - bus-pair endpoints for a bus-pair variable
# - `nothing` for a bus variable or when there is no `system_data` (flat mode).
function _graph_indices(system_data, n_dim::Int, var_name::String)
    isnothing(system_data) && return nothing
    n_dim == system_data.n_branches && return (system_data.I_branches, system_data.J_branches)
    n_dim == system_data.n_buspairs && return (system_data.I_buspairs, system_data.J_buspairs)
    n_dim == system_data.n_buses || throw(ArgumentError(
        "Variable '$var_name' has dimension $n_dim, expected $(system_data.n_branches) " *
        "(branches), $(system_data.n_buspairs) (bus pairs) or $(system_data.n_buses) (buses)"))
    return nothing
end

_add_figure_title!(fig, row::Int, title::String) =
    Label(fig[row, :], title; fontsize=20, font=:bold, halign=:center,
          tellwidth=false, tellheight=true)

"""
    viz_opf(network::Dict, variables, x, var_data...; kwargs...)
    viz_opf(system_identifier::String, variables, x, var_data...; kwargs...)

Visualize OPF solution variables and save one image per variable to `output_dir`.

Provide either a network data `Dict` similar to PowerModels format, or the name of the
system compatible with PGLib.jl, for which data is obtained with
`make_basic_network(pglib(system_identifier))`.

There are 2 calling modes, determined by the type of `var_data`:

- **Matrix mode**: `variables` is a `String`; each `var_data` is a `Matrix`
  `(n_dim × n_instances)`, one per solver.
- **Dict mode**: `variables` is a `Vector{String}`; each `var_data` is a `Dict` mapping
  variable names to matrices, one per solver.

`x` is either a single `Vector` (shared across solvers) or a `Vector` of `Vector`s (one per solver).

**Variable dispatch** (when `flat=false`):
- Dimension equals the number of **branches** → [`plot_graph_variable`](@ref), where `I`/`J` are
  `f_bus`/`t_bus` in sorted branch key order.
- Dimension equals the number of **bus pairs** → [`plot_graph_variable`](@ref), where `I`/`J` are
  `f_bus`/`t_bus` of the bus pairs in order of first occurrences in sorted branches.
- Dimension equals the number of **buses** → [`plot_variable`](@ref).

All the branches are assumed to be active and are accounted for.


**Keyword arguments**: `system_name` (label used in figure titles and output filenames; defaults to the
network's `"name"` field, or `"system"` if absent), `solver_names`, `output_dir` (default `"."`),
`vis_threshold` (default `20`), `flat` (default `false`, bypasses network loading and always uses
`plot_variable`), `xlabel` (forwarded to the underlying plotting functions; defaults to their default),
`symlog` (default `false`, draws the y-axis on a symmetric log scale),
`palette` (default `nothing`, sets the per-solver colors — a vector of colors or a `Symbol`
naming a Makie/ColorSchemes palette).
Output images are named `{system_name}_{variable}.png`.
"""
function viz_opf(
    network::Dict,
    variables::Union{String, Vector{String}},
    x,
    var_data...;
    system_name::String=get(network, "name", "system"),
    solver_names=nothing,
    output_dir::String=".",
    vis_threshold::Int=20,
    flat::Bool=false,
    xlabel=nothing,
    symlog::Bool=false,
    palette=nothing,
    alpha::Real=1.0
)
    system_data = flat ? nothing : _get_power_system_data(network)
    var_data_pairs = _var_data_pairs(variables, var_data, (2,))

    mkpath(output_dir)
    for (var_name, solvers_data) in var_data_pairs
        plot_kwargs = (; solver_names, var_name, vis_threshold, xlabel, symlog, palette, alpha)
        graph_indices = _graph_indices(system_data, size(solvers_data[1], 1), var_name)
        fig = if isnothing(graph_indices)
            plot_variable(x, solvers_data...; plot_kwargs...)
        else
            plot_graph_variable(graph_indices..., x, solvers_data...; plot_kwargs...)
        end

        _add_figure_title!(fig, 0, "$(system_name): $(var_name)")

        output_path = joinpath(output_dir, "$(system_name)_$(var_name).png")
        save(output_path, fig)
        println("Saved $output_path")
    end
end

function viz_opf(
    system_identifier::String,
    variables::Union{String, Vector{String}},
    x,
    var_data...;
    system_name::String=system_identifier,
    flat::Bool=false,
    kwargs...
)
    network = _load_network(system_identifier, flat)
    return viz_opf(network, variables, x, var_data...; system_name=system_name, flat=flat, kwargs...)
end

"""
    animate_opf(network::Dict, variables, x, time_steps, var_data...; kwargs...)
    animate_opf(system_identifier::String, variables, x, time_steps, var_data...; kwargs...)

Animate OPF solution variables over a sequence of time-stepped frames and save one
GIF per variable to `output_dir`. Animated analogue of [`viz_opf`](@ref).

Provide either a network data `Dict` similar to PowerModels format, or the name of the
system compatible with PGLib.jl, for which data is obtained with
`make_basic_network(pglib(system_identifier))`.

Two calling modes, determined by the element type of `var_data`:

- **Array mode**: `variables` is a `String`. Each `var_data` element is either an
  `(n_dim × n_instances)` `AbstractMatrix` (held constant across frames) or an
  `(n_dim × n_instances × n_frames)` 3D `AbstractArray` (animated). One element per solver.
- **Dict mode**: `variables` is a `Vector{String}`. Each `var_data` element is a
  `Dict` mapping variable names to either of the two shapes above.

Within either mode, different solvers may supply different shapes for the same variable
(e.g. an animated solver alongside a static reference) — see
[`animate_variable`](@ref) / [`animate_graph_variable`](@ref).

`x` is either a single `Vector` (shared across solvers) or a `Vector` of `Vector`s
(one per solver), following the same convention as [`viz_opf`](@ref).

`time_steps` is a `Vector` whose length must match the third dimension of any 3D `var_data`.

**Variable dispatch** (when `flat=false`):
- Dimension equals the number of **branches** → [`animate_graph_variable`](@ref),
  where `I`/`J` are `f_bus`/`t_bus` in sorted branch key order.
- Dimension equals the number of **bus pairs** → [`animate_graph_variable`](@ref),
  where `I`/`J` are `f_bus`/`t_bus` of the bus pairs in order of first occurrences
  in sorted branches.
- Dimension equals the number of **buses** → [`animate_variable`](@ref).

All the branches are assumed to be active and are accounted for.

**Keyword arguments**: `system_name` (label used in figure titles and output filenames; defaults to the
network's `"name"` field, or `"system"` if absent), `solver_names`, `output_dir` (default `"."`),
`vis_threshold` (default `20`), `flat` (default `false`, bypasses network loading and always uses
`animate_variable`), `xlabel`, `time_label` (default `"t"`), `framerate` (default `10`),
`ylims` (default `nothing`), `symlog` (default `false`, draws the y-axis on a symmetric log
scale), `palette` (default `nothing`, sets the per-solver colors — a vector of colors or a `Symbol`
naming a Makie/ColorSchemes palette).
Output files are named `{system_name}_{variable}.gif`.
"""
function animate_opf(
    network::Dict,
    variables::Union{String, Vector{String}},
    x,
    time_steps::AbstractVector,
    var_data...;
    system_name::String=get(network, "name", "system"),
    solver_names=nothing,
    output_dir::String=".",
    vis_threshold::Int=20,
    flat::Bool=false,
    xlabel=nothing,
    time_label::String="t",
    framerate::Int=10,
    ylims::Union{Nothing,Tuple{Real,Real}}=nothing,
    symlog::Bool=false,
    palette=nothing,
    alpha::Real=1.0,
)
    system_data = flat ? nothing : _get_power_system_data(network)
    var_data_pairs = _var_data_pairs(variables, var_data, (2, 3))

    mkpath(output_dir)
    n_frames = length(time_steps)
    for (var_name, solvers_data) in var_data_pairs
        plot_kwargs = (; solver_names, var_name, vis_threshold, xlabel, time_label, ylims,
                       symlog, palette, alpha)
        graph_indices = _graph_indices(system_data, size(solvers_data[1], 1), var_name)
        fig, frame_obs = if isnothing(graph_indices)
            animate_variable(x, time_steps, solvers_data...; plot_kwargs...)
        else
            animate_graph_variable(graph_indices..., x, time_steps, solvers_data...; plot_kwargs...)
        end

        _add_figure_title!(fig, -1, "$(system_name): $(var_name)")

        output_path = joinpath(output_dir, "$(system_name)_$(var_name).gif")
        record(fig, output_path, 1:n_frames; framerate=framerate) do f
            frame_obs[] = f
        end
        println("Saved $output_path")
    end
end

function animate_opf(
    system_identifier::String,
    variables::Union{String, Vector{String}},
    x,
    time_steps::AbstractVector,
    var_data...;
    system_name::String=system_identifier,
    flat::Bool=false,
    kwargs...
)
    network = _load_network(system_identifier, flat)
    return animate_opf(network, variables, x, time_steps, var_data...; system_name=system_name, flat=flat, kwargs...)
end

"""
    plot_scale_histograms_opf(system_name::String, variables, var_data...; kwargs...)
    plot_scale_histograms_opf(network::Dict, variables, var_data...; kwargs...)

Plot the per-entry scale histograms of OPF solution variables with
[`plot_scale_histograms`](@ref), one row per variable, and save the figure to `output_dir`.

The histograms do not depend on the network topology, so no network data is loaded: the first
argument only names the system, either directly or through the network's `"name"` field
(`"system"` if absent).

`variables` and `var_data` follow the two calling modes of [`viz_opf`](@ref): a `String`
with one `(n_dim × n_instances)` `Matrix` per solver, or a `Vector{String}` with one `Dict`
per solver mapping variable names to matrices. Variables may differ in dimension.

**Keyword arguments**: `system_name` (`Dict` method only; overrides the network's name),
`solver_names`, `output_dir` (default `"."`), `xlabel`, `symlog` (default `false`; a `Bool`
for every variable or a `Vector{Bool}` with one flag per variable), `bins` (default `20`),
`palette` (default `nothing`), `alpha` (default `1.0`), all forwarded to
[`plot_scale_histograms`](@ref).
The output image is named `{system_name}_scale_histograms.png`.
"""
function plot_scale_histograms_opf(
    system_name::String,
    variables::Union{String, Vector{String}},
    var_data...;
    solver_names=nothing,
    output_dir::String=".",
    xlabel=nothing,
    symlog::Union{Bool, AbstractVector{Bool}}=false,
    bins::Int=20,
    palette=nothing,
    alpha::Real=1.0
)
    var_data_pairs = _var_data_pairs(variables, var_data, (2,))
    var_names = [var_name for (var_name, _) in var_data_pairs]
    fig = plot_scale_histograms((solvers_data for (_, solvers_data) in var_data_pairs)...;
                                solver_names, var_names, xlabel, symlog, bins, palette, alpha)

    # The statistic headers occupy row 1, so row 0 places the title above them.
    _add_figure_title!(fig, 0, "$(system_name): scale histograms")

    mkpath(output_dir)
    output_path = joinpath(output_dir, "$(system_name)_scale_histograms.png")
    save(output_path, fig)
    println("Saved $output_path")
end

function plot_scale_histograms_opf(
    network::Dict,
    variables::Union{String, Vector{String}},
    var_data...;
    system_name::String=get(network, "name", "system"),
    kwargs...
)
    return plot_scale_histograms_opf(system_name, variables, var_data...; kwargs...)
end

"""
    plot_sparse_variable_opf(system_name::String, variables, x, var_data...; kwargs...)
    plot_sparse_variable_opf(network::Dict, variables, x, var_data...; kwargs...)

Plot the nonzero entries of OPF solution variables with [`plot_sparse_variable`](@ref) and save
one image per variable to `output_dir`.

The sparsity plots do not depend on the network topology, so no network data is loaded: the
first argument only names the system, either directly or through the network's `"name"` field
(`"system"` if absent).

`variables`, `x` and `var_data` follow the two calling modes of [`viz_opf`](@ref): a `String`
with one `(n_dim × n_instances)` `Matrix` per solver, or a `Vector{String}` with one `Dict`
per solver mapping variable names to matrices.

**Keyword arguments**: `system_name` (`Dict` method only; overrides the network's name),
`solver_names`, `output_dir` (default `"."`), `xlabel`, `symlog` (default `false`, draws the
color bar on a symmetric log scale), `atol` (default `1e-6`), `colormap` (default `:viridis`),
`markersize` (default `6`), `alpha` (default `1.0`), all forwarded to
[`plot_sparse_variable`](@ref).
Output images are named `{system_name}_{variable}_sparse.png`.
"""
function plot_sparse_variable_opf(
    system_name::String,
    variables::Union{String, Vector{String}},
    x,
    var_data...;
    solver_names=nothing,
    output_dir::String=".",
    xlabel=nothing,
    symlog::Bool=false,
    atol::Real=1e-6,
    colormap=:viridis,
    markersize::Real=6,
    alpha::Real=1.0
)
    var_data_pairs = _var_data_pairs(variables, var_data, (2,))

    mkpath(output_dir)
    for (var_name, solvers_data) in var_data_pairs
        fig = plot_sparse_variable(x, solvers_data...; solver_names, var_name, xlabel, symlog,
                                   atol, colormap, markersize, alpha)

        _add_figure_title!(fig, 0, "$(system_name): $(var_name)")

        output_path = joinpath(output_dir, "$(system_name)_$(var_name)_sparse.png")
        save(output_path, fig)
        println("Saved $output_path")
    end
end

function plot_sparse_variable_opf(
    network::Dict,
    variables::Union{String, Vector{String}},
    x,
    var_data...;
    system_name::String=get(network, "name", "system"),
    kwargs...
)
    return plot_sparse_variable_opf(system_name, variables, x, var_data...; kwargs...)
end
