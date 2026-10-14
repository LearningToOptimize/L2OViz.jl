# Per-entry values across all instances (and frames for 3D animation data).
entry_values(d::AbstractMatrix, k::Int) = @view d[k, :]
entry_values(d::AbstractArray{<:Any, 3}, k::Int) = vec(@view d[k, :, :])

# Return the y-data to pass to `scatter!` for a single (entry, solver) pair.
# A Matrix becomes a plain Vector (does not relate to any animation frame).
# A 3D array becomes an Observable that re-fires whenever `frame_obs` changes.
plot_y(d::AbstractMatrix, k::Int, ::Any) = d[k, :]
plot_y(d::AbstractArray{<:Any, 3}, k::Int, frame_obs::Observable) = @lift d[k, :, $frame_obs]

# Resolve a set of names: default to `"$prefix $i"` placeholders if not supplied, otherwise
# validate length. `kind` names the keyword argument in error messages (e.g. "solver_names").
resolve_names(::Nothing, n::Int, ::AbstractString, prefix::AbstractString) =
    ["$prefix $i" for i in 1:n]
function resolve_names(names, n::Int, kind::AbstractString, ::AbstractString)
    if isa(names, AbstractVector{<:AbstractString})
        length(names) == n || throw(ArgumentError(
            "$kind must have length $n, got $(length(names))"))
        return collect(names)
    else
        throw(ArgumentError("$kind must be a Vector of Strings or nothing"))
    end
end

# Resolve `solver_names`: default to placeholder names if not supplied, otherwise validate length.
resolve_solver_names(solver_names, n_solvers::Int) =
    resolve_names(solver_names, n_solvers, "solver_names", "Solver")

function validate_var_data_dims(var_data)
    # All matrices must have the same number of rows (dimension of the variable)
    n_entries = size(var_data[1], 1)
    n_entries >= 1 ||
        throw(ArgumentError("var_data must have at least one entry (size(var_data[1], 1) > 0)"))
    for (i, d) in enumerate(var_data)
        size(d, 1) == n_entries || throw(DimensionMismatch(
            "Variable dimension of solver $(i): $(size(d, 1)); mismatches with variable dimension of solver 1: $(n_entries)"))
    end
    return n_entries
end

function validate_var_data_dims(var_data, n_frames::Int)
    # All arrays must have the same length in the first dimension (dimension of the variable)
    # 3D arrays must match the number of frames.
    # (Matrices have no time dimension and are constant across frames.)
    n_entries = validate_var_data_dims(var_data)
    for (i, d) in enumerate(var_data)
        if ndims(d) == 3
            size(d, 3) == n_frames || throw(DimensionMismatch(
                "Number of frames in data of solver $(i): $(size(d, 3)); does not equal length(time_steps) = $(n_frames)"))
        end
    end
    return n_entries
end

# Resolve `x` into one Vector per solver. `x` may be a single shared Vector (broadcast to
# every solver) or a Vector of per-solver Vectors. Validates that each `x` Vector's length
# matches the number of instances in that solver's data array.
function resolve_x_vecs(x, var_data)
    n_solvers = length(var_data)
    if all(isa.(x, AbstractVector))
        # x is multiple vectors, one per solver
        x_vecs = collect(x)
        length(x_vecs) == n_solvers || throw(ArgumentError(
            "Number of x vectors ($(length(x_vecs))) must equal number of data arrays ($n_solvers)"))
        for (i, d) in enumerate(var_data)
            n_instances_i = size(d, 2)
            length(x_vecs[i]) == n_instances_i || throw(DimensionMismatch(
                "Length of x[$(i)] ($(length(x_vecs[i]))) must equal number of columns in data array $(i) ($n_instances_i)"))
        end
        return x_vecs
    elseif isa(x, AbstractVector)
        # x is a single vector shared across all solvers
        for (i, d) in enumerate(var_data)
            n_instances_i = size(d, 2)
            length(x) == n_instances_i || throw(DimensionMismatch(
                "Length of x ($(length(x))) must equal number of columns in data array $(i) ($n_instances_i)"))
        end
        return [x for _ in 1:n_solvers]
    else
        throw(ArgumentError("x must be either a single AbstractVector or multiple AbstractVectors (one per solver); got $(typeof(x))"))
    end
end

# Per-solver scatter color, cycling through `palette` if there are more solvers than palette
# entries. `palette` may be:
# - `nothing`: Makie's default categorical palette (wong_colors).
# - a `Symbol`: the name of a Makie/ColorSchemes palette (e.g. `:tab10`, `:viridis`).
# - a vector of colors (anything Makie accepts, e.g. `[:red, :blue]`).
function solver_palette(n_solvers::Int, palette=nothing)
    colors = resolve_palette_colors(palette, n_solvers)
    isempty(colors) && throw(ArgumentError("palette must contain at least one color"))
    return [colors[mod1(i, length(colors))] for i in 1:n_solvers]
end

resolve_palette_colors(::Nothing, ::Int) = Makie.wong_colors()
resolve_palette_colors(palette, ::Int) = palette

# Resolve a named palette into a vector of colors. `categorical_colors` samples continuous
# colormaps (e.g. `:viridis`) into `n_solvers` evenly spaced colors so distinct solvers stay
# visually separated, and returns categorical schemes' (e.g. `:tab10`) discrete swatches. It
# errors only when a categorical scheme has fewer swatches than `n_solvers`, in which case we
# fall back to the full swatch set for `solver_palette` to cycle through.
function resolve_palette_colors(palette::Symbol, n_solvers::Int)
    try
        return Makie.categorical_colors(palette, n_solvers)
    catch err
        err isa ErrorException || rethrow()
        return Makie.to_colormap(palette)
    end
end

# Per-row significance score for the variable: aggregate every solver's values at row `k`
# across all instances (and all frames, for 3D data) and reduce with `significance_fn`.
# Works for both static and animated input transparently via `entry_values` dispatch.
function compute_entry_scores(var_data, n_entries::Int, significance_fn)
    return [significance_fn(vcat([Vector(entry_values(d, k)) for d in var_data]...))
            for k in 1:n_entries]
end

# Min/max over the given entries of every solver's data (and over all frames for 3D data).
function data_extrema(var_data, entry_indices)
    ymin = Inf
    ymax = -Inf
    for k in entry_indices
        for d in var_data
            slice = entry_values(d, k)
            ymin = min(ymin, minimum(slice))
            ymax = max(ymax, maximum(slice))
        end
    end
    return ymin, ymax
end

# Auto y-axis limits with 5% padding, computed across the given entries of every solver's
# data. Used by animate_* to fix a single y-range across all frames so that the animation
# does not auto-rescale frame-by-frame.
# The padding is applied in the axis' transformed space, so that a nonlinear `yscale` gets
# visually even margins.
function compute_ylim_range(var_data, entry_indices, yscale)
    ymin, ymax = data_extrema(var_data, entry_indices)
    lo, hi = yscale(ymin), yscale(ymax)
    span = hi - lo
    pad = span > 0 ? 0.05 * span : 0.5 * max(abs(hi), 1.0)
    inverse = Makie.inverse_transform(yscale)
    return inverse(lo - pad), inverse(hi + pad)
end

# Apply user-supplied ylims, or compute & apply auto ylims with padding.
# Since the axes are linked, setting limits on one axis propagates to every panel.
function apply_fixed_ylims!(ax, ylims_user, var_data, entry_indices, yscale)
    if isnothing(ylims_user)
        lo, hi = compute_ylim_range(var_data, entry_indices, yscale)
        ylims!(ax, lo, hi)
    else
        ylims!(ax, ylims_user[1], ylims_user[2])
    end
end

# Largest number of decades a symlog axis is allowed to span on each side of zero. Bounds the
# linear region from below so that a single near-zero value cannot stretch the axis out.
const SYMLOG_MAX_DECADES = 8

# Half-width of the linear region of a symlog axis, derived from a flat iterable of values:
# the smallest nonzero magnitude present, floored at `SYMLOG_MAX_DECADES` below the largest
# magnitude. Returns `nothing` when every value is zero, where a symlog axis is not meaningful.
function symlog_linthresh_from_values(values)
    max_abs = 0.0
    min_nonzero_abs = Inf
    for value in values
        magnitude = abs(value)
        magnitude == 0 && continue
        max_abs = max(max_abs, magnitude)
        min_nonzero_abs = min(min_nonzero_abs, magnitude)
    end
    max_abs > 0 || return nothing
    return clamp(min_nonzero_abs, max_abs / 10.0^SYMLOG_MAX_DECADES, max_abs)
end

# Axis scale for the given entries of every solver's data (and all frames for 3D data).
resolve_yscale(symlog::Bool, var_data, entry_indices) = resolve_symlog_transform(
    symlog, (value for k in entry_indices for d in var_data for value in entry_values(d, k)))

# Resolve a symlog axis transform from a flat iterable of values, falling back to `identity`
# when a symlog scale is not meaningful (every value is zero).
function resolve_symlog_transform(symlog::Bool, values)
    symlog || return identity
    linthresh = symlog_linthresh_from_values(values)
    if isnothing(linthresh)
        @warn "symlog=true, but every value is zero; falling back to a linear axis."
        return identity
    end
    return Makie.Symlog10(linthresh)
end

# Histogram bin edges (length `n_bins + 1`) spanning the extrema of `values`. For a nonlinear
# `transform` the edges are evenly spaced in transformed space and inverted back to data
# space, so that an axis using `xscale=transform` renders bars at even visual width. A
# degenerate range (all values equal) is padded symmetrically so the edges stay increasing.
function scale_histogram_bins(values, n_bins::Int, transform)
    lo, hi = extrema(values)
    tlo, thi = transform(lo), transform(hi)
    if !(thi > tlo)
        pad = tlo == 0 ? 0.5 : 0.5 * abs(tlo)
        tlo -= pad
        thi += pad
    end
    inverse = Makie.inverse_transform(transform)
    edges = [inverse(tlo + (thi - tlo) * i / n_bins) for i in 0:n_bins]
    edges[end] = max(edges[end], nextfloat(hi))
    return edges
end

# Resolve the `symlog` argument of `plot_scale_histograms` into one flag per variable. A single
# `Bool` applies to every variable; a `Vector{Bool}` selects the scale per variable.
resolve_symlog_flags(symlog::Bool, n_vars::Int) = fill(symlog, n_vars)
function resolve_symlog_flags(symlog::AbstractVector, n_vars::Int)
    all(x -> x isa Bool, symlog) ||
        throw(ArgumentError("symlog must be a Bool or a Vector of Bools"))
    length(symlog) == n_vars || throw(ArgumentError(
        "symlog must have length $n_vars (one per variable), got $(length(symlog))"))
    return collect(Bool, symlog)
end
resolve_symlog_flags(symlog, ::Int) =
    throw(ArgumentError("symlog must be a Bool or a Vector of Bools, got $(typeof(symlog))"))

# Panel titles. Vector variant produces "[k]" or "var_name[k]"; graph variant produces
# "[i,j]" or "var_name[i,j]".
entry_label(var_name::AbstractString, idx::Int) =
    isempty(var_name) ? "[$idx]" : "$(var_name)[$idx]"
entry_label(var_name::AbstractString, i::Int, j::Int) =
    isempty(var_name) ? "[$i,$j]" : "$(var_name)[$i,$j]"

# Vector layout: arrange `n_plot` panels into a near-square grid (n_cols ≈ sqrt(n_plot)).
function vector_grid_layout(n_plot::Int)
    n_cols = ceil(Int, sqrt(n_plot))
    n_rows = ceil(Int, n_plot / n_cols)
    return n_rows, n_cols
end

# Grid cell `(row, col)` of the `k`-th panel in a row-major vector layout with `n_cols` columns.
vector_grid_cell(k::Int, n_cols::Int) = (div(k - 1, n_cols) + 1, mod(k - 1, n_cols) + 1)

# Compute dimensions `(n_rows, n_cols)` of the plot grid and the mapping (i, j) → grid cell.
# Rows and columns share one sorted vertex set (square grid) unless `bipartite`, in which case
# rows are indexed by the unique `I_plot` and columns by the unique `J_plot`.
function matrix_grid_layout(I_plot, J_plot; bipartite::Bool=false)
    # (i, j) should be sorted in the resulting grid
    if bipartite
        row_vertices_sorted = sort(unique(I_plot))
        col_vertices_sorted = sort(unique(J_plot))
    else
        row_vertices_sorted = col_vertices_sorted = sort(union(unique(I_plot), unique(J_plot)))
    end
    row_index_map = Dict(v => idx for (idx, v) in enumerate(row_vertices_sorted))
    col_index_map = Dict(v => idx for (idx, v) in enumerate(col_vertices_sorted))
    return length(row_vertices_sorted), length(col_vertices_sorted),
           (i, j) -> (row_index_map[i], col_index_map[j])
end

# Draw the vector-variable panel grid into `parent` (a Figure or GridLayout). Returns
# (axes, legend_handles). `frame_obs` is `nothing` for static plots and an Observable{Int}
# for animations; this is what selects between a plain Vector and a lifted Observable in
# `plot_y` for each 3D solver array.
function draw_vector_panels!(parent, var_data, x_vecs, x_label, var_name,
                             selected_indices, solver_colors, n_cols; frame_obs=nothing,
                             yscale=identity, alpha=1.0)
    legend_handles = []
    axes = Axis[]
    for (k, data_idx) in enumerate(selected_indices)
        grid_row, grid_col = vector_grid_cell(k, n_cols)
        ax = Axis(parent[grid_row, grid_col];
                  title=entry_label(var_name, data_idx), xlabel=x_label, yscale=yscale)
        push!(axes, ax)
        for (i, d) in enumerate(var_data)
            p = scatter!(ax, x_vecs[i], plot_y(d, data_idx, frame_obs);
                         color=solver_colors[i], alpha=alpha)
            k == 1 && push!(legend_handles, p)
        end
    end
    return axes, legend_handles
end

# Same as `draw_vector_panels!`, but for a graph variable. `grid_pos` maps the original (I, J)
# to (row, col) cells of `gl` (a GridLayout). `selected_indices` maps the k-th plotted entry
# back to its row in each solver's data array.
function draw_matrix_panels!(gl, var_data, x_vecs, x_label, var_name,
                             I_plot, J_plot, selected_indices, grid_pos,
                             solver_colors; frame_obs=nothing, yscale=identity, alpha=1.0)
    legend_handles = []
    axes = Axis[]
    for k in eachindex(I_plot)
        gr, gc = grid_pos(I_plot[k], J_plot[k])
        ax = Axis(gl[gr, gc];
                  title=entry_label(var_name, I_plot[k], J_plot[k]), xlabel=x_label,
                  yscale=yscale)
        push!(axes, ax)
        coo_row = selected_indices[k]
        for (i, d) in enumerate(var_data)
            p = scatter!(ax, x_vecs[i], plot_y(d, coo_row, frame_obs);
                         color=solver_colors[i], alpha=alpha)
            k == 1 && push!(legend_handles, p)
        end
    end
    return axes, legend_handles
end

# Each solver's tuple of (real, imag) matrices must have matching sizes.
function validate_polar_pair_dims(var_data)
    for (i, (real_part, imag_part)) in enumerate(var_data)
        size(real_part) == size(imag_part) || throw(DimensionMismatch(
            "Real and imaginary parts of solver $(i) have different sizes: $(size(real_part)) vs $(size(imag_part))"))
    end
end

# Convert Cartesian (real, imag) components to polar (radius, angle in radians).
polar_coordinates(real_part::AbstractMatrix, imag_part::AbstractMatrix) =
    hypot.(real_part, imag_part), atan.(imag_part, real_part)

# Marker size range used to encode the instance parameter `x` in polar plots.
const POLAR_MIN_MARKERSIZE = 4.0
const POLAR_MAX_MARKERSIZE = 18.0

# Map each solver's `x` linearly onto [POLAR_MIN_MARKERSIZE, POLAR_MAX_MARKERSIZE] using the
# global x range across solvers, so sizes are comparable between solvers. Returns the per-solver
# size vectors and `(x value, marker size)` pairs for the legend key.
function x_marker_sizes(x_vecs; n_legend_entries::Int=3)
    x_min = minimum(minimum, x_vecs)
    x_max = maximum(maximum, x_vecs)
    size_span = POLAR_MAX_MARKERSIZE - POLAR_MIN_MARKERSIZE
    if x_max > x_min
        to_size = v -> POLAR_MIN_MARKERSIZE + size_span * (v - x_min) / (x_max - x_min)
        legend_x_values = range(x_min, x_max, length=n_legend_entries)
    else
        to_size = _ -> POLAR_MIN_MARKERSIZE + size_span / 2
        legend_x_values = [x_min]
    end
    marker_sizes = [to_size.(x_vec) for x_vec in x_vecs]
    return marker_sizes, [(v, to_size(v)) for v in legend_x_values]
end

# Legend below the polar panels: one group for solver colors and one for the marker-size key
# of the instance parameter.
function add_polar_legend!(legend_cell, solver_colors, solver_names, size_key, x_label)
    solver_handles = [MarkerElement(marker=:circle, color=color,
                                    markersize=(POLAR_MIN_MARKERSIZE + POLAR_MAX_MARKERSIZE) / 2)
                      for color in solver_colors]
    size_handles = [MarkerElement(marker=:circle, color=:gray50, markersize=marker_size)
                    for (_, marker_size) in size_key]
    size_labels = [string(round(x_value; sigdigits=3)) for (x_value, _) in size_key]
    Legend(legend_cell, [solver_handles, size_handles], [solver_names, size_labels],
           ["Solver", x_label]; orientation=:horizontal, tellwidth=false, tellheight=true)
end

# Radial ticks for a PolarAxis whose radii were transformed by `rscale`: ticks are placed in
# transformed space but labeled with untransformed radii. Labels are always `RichText`, because
# PolarAxis fixes the label type at construction and errors if the automatic formatter later
# switches between `String` and `RichText` (e.g. scientific notation for large radii).
function polar_rticks(rscale)
    inverse = Makie.inverse_transform(rscale)
    return function (rmin_scaled, rmax_scaled)
        tick_radii, tick_labels = polar_tick_radii_and_labels(rscale, inverse(rmin_scaled),
                                                              inverse(rmax_scaled))
        return rscale.(tick_radii),
               [label isa Makie.RichText ? label : rich(label) for label in tick_labels]
    end
end

polar_tick_radii_and_labels(rscale, rmin, rmax) =
    Makie.get_ticks(Makie.automatic, rscale, Makie.automatic, rmin, rmax)

# Maximum number of ticks (besides 0) on a symlog radial axis.
const POLAR_MAX_SYMLOG_TICKS = 4

# Symlog radii: 0 plus log-spaced ticks between the linear threshold and `rmax`, thinned by a
# uniform stride to at most POLAR_MAX_SYMLOG_TICKS. Decade ticks (10^k) are used when at least
# two decades fit; otherwise 1-2-5 ticks fill the narrower range. Unlike Makie's `SymlogTicks`,
# decades below 1 are allowed, so small radii also get log-spaced ticks.
function polar_tick_radii_and_labels(rscale::Makie.Symlog10, rmin, rmax)
    linthresh = rscale.upper
    exponents = floor(Int, log10(linthresh)):floor(Int, log10(rmax))
    decade_radii = filter(radius -> linthresh <= radius <= rmax, exp10.(exponents))
    use_decades = length(decade_radii) >= 2
    candidate_radii = use_decades ? decade_radii :
        filter(radius -> linthresh <= radius <= rmax,
               [mantissa * exp10(k) for k in exponents for mantissa in (1, 2, 5)])
    length(candidate_radii) >= 2 || return Makie.get_ticks(Makie.automatic, identity,
                                                           Makie.automatic, rmin, rmax)

    # Keep the largest candidate and step down from it, so the outermost tick is closest to the rim.
    stride = cld(length(candidate_radii), POLAR_MAX_SYMLOG_TICKS)
    shown_radii = reverse(reverse(candidate_radii)[1:stride:end])
    shown_labels = use_decades ?
        [rich("10", superscript(string(round(Int, log10(radius))), offset=Vec2f(0.1, 0.0)))
         for radius in shown_radii] :
        Makie.get_ticklabels(Makie.automatic, shown_radii)
    return [0.0; shown_radii], [rich("0"); shown_labels]
end

# Draw one PolarAxis per panel into `parent` (a Figure or GridLayout). Panel `k` is placed at
# `panel_cells[k]` and shows row `data_indices[k]` of every solver's data. Radii are plotted after
# applying `rscale` (PolarAxis has no radial scale attribute). Every panel uses the transformed
# radius range [0, rmax_scaled], since PolarAxis cannot be linked.
function draw_polar_panels!(parent, radius_data, angle_data, marker_sizes, panel_cells,
                            panel_titles, data_indices, solver_colors, rmax_scaled;
                            rscale=identity, alpha=1.0)
    axes = PolarAxis[]
    for (k, data_idx) in enumerate(data_indices)
        grid_row, grid_col = panel_cells[k]
        ax = PolarAxis(parent[grid_row, grid_col]; title=panel_titles[k],
                       rticks=polar_rticks(rscale), rtickangle=π / 2)
        push!(axes, ax)
        for (i, (radius, angle)) in enumerate(zip(radius_data, angle_data))
            scatter!(ax, angle[data_idx, :], rscale.(radius[data_idx, :]);
                     markersize=marker_sizes[i], color=solver_colors[i], alpha=alpha)
        end
        rlims!(ax, 0, rmax_scaled)
    end
    return axes
end
