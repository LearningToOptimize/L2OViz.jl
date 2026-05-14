"""Default significance: 1-norm (sum of absolute values) across problem instances."""
default_significance(values::AbstractVector) = sum(abs, values)

"""
    select_variable_entries(scores, vis_threshold) -> (selected_indices, filtered::Bool)

Return the indices of the `vis_threshold` highest-scoring entries, sorted ascending.
If `length(scores) ≤ vis_threshold`, returns all indices and `filtered = false`.
"""
function select_variable_entries(scores::Vector{<:Real}, vis_threshold::Int)
    n = length(scores)
    n <= vis_threshold && return collect(1:n), false
    return sort(partialsortperm(scores, 1:vis_threshold, rev=true)), true
end

"""
    select_matrix_entries(entry_scores, I, J, vis_threshold, symmetric)
        -> (I_plot, J_plot, data_row_indices, sel_rows, sel_cols, filtered::Bool)

Determine which matrix entries to display, with optional thresholding and symmetric mirroring.

Row/column significance = max entry score among entries in that row/column.
Filtering selects the top-vis_threshold rows and columns (independently when `symmetric = false`).

When `symmetric = true`, off-diagonal entries in the (filtered) COO are mirrored so the full
symmetric subgrid is shown. `data_row_indices[k]` maps back into the original COO rows of
var_data; mirrored entries reuse the same index as their originals.
"""
function select_matrix_entries(entry_scores::Vector{<:Real},
                                I::Vector{Int}, J::Vector{Int},
                                vis_threshold::Union{Int, Nothing},
                                symmetric::Bool)
    n_nonzeros = length(I)
    needs_filter = !isnothing(vis_threshold) && n_nonzeros > vis_threshold

    if symmetric
        # Build per-index significance: max score of all entries where the index appears
        all_indices = sort(unique(vcat(I, J)))
        index_scores = Dict{Int, Float64}()
        for k in 1:n_nonzeros
            s = entry_scores[k]
            index_scores[I[k]] = max(get(index_scores, I[k], 0.0), s)
            index_scores[J[k]] = max(get(index_scores, J[k], 0.0), s)
        end

        if needs_filter
            k_sel = floor(Int, vis_threshold)
            sorted_idx = sort(all_indices, by=i -> index_scores[i], rev=true)
            sel_indices = sort(sorted_idx[1:min(k_sel, end)])
        else
            sel_indices = all_indices
        end

        sel_set = Set(sel_indices)
        orig_idx = findall(k -> I[k] ∈ sel_set && J[k] ∈ sel_set, 1:n_nonzeros)

        return I[orig_idx], J[orig_idx], orig_idx, sel_indices, sel_indices, needs_filter

    else
        # Score rows and columns independently using a single O(n_nonzeros) pass
        row_scores = Dict{Int, Float64}()
        col_scores = Dict{Int, Float64}()
        for k in 1:n_nonzeros
            s = entry_scores[k]
            row_scores[I[k]] = max(get(row_scores, I[k], 0.0), s)
            col_scores[J[k]] = max(get(col_scores, J[k], 0.0), s)
        end

        unique_rows = sort(collect(keys(row_scores)))
        unique_cols = sort(collect(keys(col_scores)))

        if !needs_filter
            return I, J, collect(1:n_nonzeros), unique_rows, unique_cols, false
        end

        k_sel = floor(Int, vis_threshold)
        sel_rows = sort(sort(unique_rows, by=r -> row_scores[r], rev=true)[1:min(k_sel, end)])
        sel_cols = sort(sort(unique_cols, by=c -> col_scores[c], rev=true)[1:min(k_sel, end)])
        sel_row_set, sel_col_set = Set(sel_rows), Set(sel_cols)

        orig_idx = findall(k -> I[k] ∈ sel_row_set && J[k] ∈ sel_col_set, 1:n_nonzeros)
        return I[orig_idx], J[orig_idx], orig_idx, sel_rows, sel_cols, true
    end
end
