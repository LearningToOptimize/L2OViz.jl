# exp/viz_opf.jl
# Utility script for visualizing OPF solution data from PowerModels.jl / PGLib.jl.
#
# Usage:
#   julia --project=exp exp/viz_opf.jl \
#       --system pglib_opf_case14_ieee \
#       --h5 solver_a.h5 solver_b.h5 \
#       --solver-names "Solver A" "Solver B" \
#       --variables pg qg \
#       --output-dir /path/to/output \
#       [--flat]
#
# Required packages in the project environment: ArgParse, HDF5, L2OViz, PGLib, PowerModels.

using ArgParse
using CairoMakie
using HDF5
using L2OViz
using PGLib
using PowerModels

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--system"
            help = "PGLib system name (e.g. pglib_opf_case14_ieee)"
            required = true
        "--h5"
            help = "Path to the .h5 data file for each solver (one or more)"
            nargs = '+'
            required = true
        "--solver-names"
            help = "Display name for each solver, in the same order as --h5 (defaults to 'Solver 1', 'Solver 2', ...)"
            nargs = '+'
        "--variables"
            help = "Names of variables to visualize; one image is saved per variable"
            nargs = '+'
            required = true
        "--output-dir"
            help = "Directory where output images are saved (default: current directory)"
            default = "."
        "--flat"
            help = "Use plot_variable for all variables instead of plot_matrix_variable"
            action = :store_true
        "--vis-threshold"
            help = "Maximum number of entries (plot_variable) or rows/columns (plot_matrix_variable) to visualize per variable"
            arg_type = Int
            default = 20
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()

    system_name   = args["system"]
    h5_paths      = args["h5"]
    variables     = args["variables"]
    output_dir    = args["output-dir"]
    use_flat_plot  = args["flat"]
    vis_threshold  = args["vis-threshold"]
    solver_names  = isempty(args["solver-names"]) ?
                        ["Solver $i" for i in 1:length(h5_paths)] :
                        args["solver-names"]
    @assert length(solver_names) == length(h5_paths) (
        "Number of --solver-names ($(length(solver_names))) must match " *
        "number of --h5 paths ($(length(h5_paths)))"
    )

    # Load network to get branch COO (f_bus → I, t_bus → J, sorted by branch key)
    # and bus/branch counts for dimension-based dispatch. Skipped in --flat mode.
    if !use_flat_plot
        network = make_basic_network(pglib(system_name))
        branch_dict = network["branch"]
        sorted_branch_keys = sort(collect(keys(branch_dict)), by=k -> parse(Int, string(k)))
        I_branches = [branch_dict[k]["f_bus"] for k in sorted_branch_keys]
        J_branches = [branch_dict[k]["t_bus"] for k in sorted_branch_keys]
        n_branches = length(sorted_branch_keys)
        n_buses    = length(network["bus"])
    end

    # Open all solver h5 files; read datasets directly without materializing full Dicts.
    # Each file must contain one dataset per variable, shaped (n_dim, n_instances),
    # and "pd" and "qd" datasets shaped (n_buses, n_instances).
    h5_files = [h5open(path, "r") for path in h5_paths]
    try
        # x-axis per solver: total load (active + reactive demand summed across buses)
        # per instance. Each solver may cover different problem instances.
        x_per_solver = [vec(sum(read(f["pd"]) .+ read(f["qd"]), dims=1)) for f in h5_files]

        mkpath(output_dir)
        for var_name in variables
            var_matrices = [read(f[var_name]) for f in h5_files]

            n_dim = size(var_matrices[1], 1)

            if use_flat_plot
                fig = plot_variable(x_per_solver, var_matrices...;
                                    solver_names=solver_names,
                                    var_name=var_name,
                                    xlabel="Total Load (pd + qd)",
                                    vis_threshold=vis_threshold)
            elseif n_dim == n_branches
                fig = plot_matrix_variable(I_branches, J_branches, x_per_solver, var_matrices...;
                                           solver_names=solver_names,
                                           var_name=var_name,
                                           xlabel="Total Load (pd + qd)",
                                           vis_threshold=vis_threshold)
            else
                @assert n_dim == n_buses (
                    "Variable '$var_name' has dimension $n_dim, expected $n_branches (branches) or $n_buses (buses)"
                )
                fig = plot_variable(x_per_solver, var_matrices...;
                                    solver_names=solver_names,
                                    var_name=var_name,
                                    xlabel="Total Load (pd + qd)",
                                    vis_threshold=vis_threshold)
            end

            output_path = joinpath(output_dir, "$(system_name)_$(var_name).png")
            save(output_path, fig)
            println("Saved $output_path")
        end
    finally
        foreach(close, h5_files)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
