"""Default significance: 1-norm (sum of absolute values) across problem instances."""
default_significance(values::AbstractVector) = sum(abs, values)

"""
    select_variable_entries(scores, vis_threshold) -> (selected_indices, filtered::Bool)

Return the indices of the `vis_threshold` highest-scoring entries, sorted ascending.
If `length(scores) ≤ vis_threshold`, returns all indices and `filtered = false`.
"""
function select_variable_entries(scores::Vector{<:Real}, vis_threshold::Int)
    vis_threshold > 0 || throw(ArgumentError("vis_threshold must be positive"))
    length(scores) <= vis_threshold && return collect(1:length(scores)), false
    return sort(partialsortperm(scores, 1:vis_threshold, rev=true)), true
end

"""
    dedup_repeated_entries(I_plot, J_plot, selected_indices, entry_scores)
        -> (I_plot, J_plot, selected_indices)

Collapse plotted entries that share the exact same `(i, j)` coordinate, keeping only the
occurrence with the highest `entry_scores` value and emitting a warning for each colliding
coordinate. This handles variables whose COO is a multi-edge graph rather than a true matrix.

`selected_indices[k]` is the original COO row of the `k`-th plotted entry, so the tie-break
compares `entry_scores[selected_indices[k]]`. The surviving entries' original rows are returned
in the (filtered) `selected_indices`.
"""
function dedup_repeated_entries(I_plot::Vector{Int}, J_plot::Vector{Int},
                                selected_indices::Vector{Int}, entry_scores::Vector{<:Real})
    # For each ordered coordinate pair, remember the index of its highest-scoring occurrence so far
    best_idx_for_pair = Dict{Tuple{Int, Int}, Int}()
    # Coordinates seen more than once; warned about after the scan.
    duplicate_pairs = Set{Tuple{Int, Int}}()
    for k in eachindex(I_plot)
        pair = (I_plot[k], J_plot[k])
        previous_position = get(best_idx_for_pair, pair, nothing)
        if isnothing(previous_position)
            # no duplicate seen so far
            best_idx_for_pair[pair] = k
        else
            push!(duplicate_pairs, pair)
            # keep the higher-scoring occurrence
            if entry_scores[selected_indices[k]] > entry_scores[selected_indices[previous_position]]
                best_idx_for_pair[pair] = k
            end
        end
    end

    isempty(duplicate_pairs) && return I_plot, J_plot, selected_indices

    for pair in collect(duplicate_pairs)
        @warn "$pair is repeated in the selected coordinates; keeping only the highest-scoring occurrence. Suggest using flat visualization."
    end

    # Keep one index per pair, preserving the original relative ordering of the entries.
    kept_indices = sort(collect(values(best_idx_for_pair)))
    return I_plot[kept_indices], J_plot[kept_indices], selected_indices[kept_indices]
end

"""
    select_matrix_entries(entry_scores, I, J, vis_threshold)
        -> (I_plot, J_plot, selected_indices, filtered::Bool)

Determine which entries of a symmetric matrix to display, with optional thresholding.

The matrix is assumed symmetric and the COO coordinates should not contain both `(i, j)` and
`(j, i)`. If number of unique rows/columns exceeds `vis_threshold`, select a principal submatrix
on the `vis_threshold` highest-scoring rows/columns.

Returns the coordinates `I_plot`, `J_plot` of the entries to display, together with
`selected_indices`, their positions in the original `I`/`J` (so `I_plot == I[selected_indices]`
and each entry's data is row `selected_indices[k]` of every solver's `var_data`). Both are
needed downstream: `(I_plot[k], J_plot[k])` places panel `k` on the grid, while
`selected_indices[k]` locates its data. `filtered` indicates whether thresholding was applied.

If the entries selected for display contain repeated exact `(i, j)` coordinates (which can
happen when the variable encodes a multi-edge graph), only the highest-scoring occurrence of
each coordinate is kept and a warning is emitted per colliding coordinate (see
`dedup_repeated_entries`).
"""
function select_matrix_entries(entry_scores::Vector{<:Real},
                                I::Vector{Int}, J::Vector{Int},
                                vis_threshold::Int)
    vis_threshold > 0 || throw(ArgumentError("vis_threshold must be positive"))

    # Since the COO is a half representation, both row and column index sets of the full
    # symmetric matrix are union(unique(I), unique(J)).
    nz_col_indices = sort(union(unique(I), unique(J)))

    # Early termination: no need to filter (but still collapse any repeated coordinates).
    if length(nz_col_indices) <= vis_threshold
        I_plot, J_plot, selected_indices =
            dedup_repeated_entries(I, J, collect(1:length(I)), entry_scores)
        return I_plot, J_plot, selected_indices, false
    end

    # Score each index by the maximum of the scores of its entries. Since we assume I and J
    # do not contain both (i, j) and (j, i), we need to aggregate the scores over both rows
    # and columns of the input matrix T, to effectively aggregate over columns of the full
    # symmetric matrix.
    col_scores = Dict{Int, Float64}()
    for k in 1:length(I)
        s = entry_scores[k]
        col_scores[I[k]] = max(get(col_scores, I[k], -Inf), s)  # max over the rows of T
        col_scores[J[k]] = max(get(col_scores, J[k], -Inf), s)  # max over the columns of T
    end
    # Select the vis_threshold highest-scoring rows/columns of the full symmetric matrix.
    sorted_indices = sort(nz_col_indices, by=c -> col_scores[c], rev=true)
    selected_vertices = Set(sorted_indices[1:min(vis_threshold, end)])
    # Keep the entries whose endpoints are both among the selected rows/columns. These indices
    # point into the original I/J (and thus into each solver's data rows).
    selected_indices = findall(
        k -> I[k] ∈ selected_vertices && J[k] ∈ selected_vertices,
        1:length(I)
    )
    isempty(selected_indices) && throw(ArgumentError("No entries selected for display; consider increasing vis_threshold"))
    # Collapse any duplicate (i, j) coordinates among the selected entries.
    I_plot, J_plot, selected_indices = dedup_repeated_entries(
        I[selected_indices], J[selected_indices],
        selected_indices, entry_scores
    )
    return I_plot, J_plot, selected_indices
end
