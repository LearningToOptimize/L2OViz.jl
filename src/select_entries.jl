"""Default significance: 1-norm (sum of absolute values) across problem instances."""
default_significance(values::AbstractVector) = sum(abs, values)

"""
    select_variable_entries(scores, vis_threshold) -> selected_indices

Return the indices of the `vis_threshold` highest-scoring entries, sorted ascending.
If `length(scores) ≤ vis_threshold`, returns all indices.
"""
function select_variable_entries(scores::Vector{<:Real}, vis_threshold::Int)
    vis_threshold > 0 || throw(ArgumentError("vis_threshold must be positive"))
    length(scores) <= vis_threshold && return collect(1:length(scores))
    return sort(partialsortperm(scores, 1:vis_threshold, rev=true))
end

"""
    dedup_edges_by_scores(I_plot, J_plot, selected_indices, entry_scores)
        -> (I_plot, J_plot, selected_indices)

Collapse repeated vertex pair `(i, j)` in `I_plot`, `J_plot`, keeping only the occurrence
with the highest `entry_scores` value and emitting a warning for each colliding coordinate.
"""
function dedup_edges_by_scores(I_plot::Vector{Int}, J_plot::Vector{Int},
                                selected_indices::Vector{Int}, entry_scores::Vector{<:Real})
    # For each ordered coordinate pair, remember the index of its highest-scoring occurrence so far
    best_idx_for_pair = Dict{Tuple{Int, Int}, Int}()
    # Coordinates seen more than once; warned about after the scan.
    duplicate_pairs = Set{Tuple{Int, Int}}()
    for k in eachindex(I_plot)
        pair = (I_plot[k], J_plot[k])
        previous_position = get(best_idx_for_pair, pair, 0)
        if previous_position == 0
            # no duplicate seen so far
            best_idx_for_pair[pair] = k
        else
            push!(duplicate_pairs, pair)
            # keep the higher-scoring occurrence
            if entry_scores[k] > entry_scores[previous_position]
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

# Select the `vis_threshold` vertices in `V` with the highest scores, where each vertex is
# scored by the maximum score of its incident edges (for each vector `endpoints` in
# `endpoint_vectors`, `endpoints[k]` is an endpoint of the edge with score `entry_scores[k]`).
function top_scoring_vertices(V::Vector{Int}, endpoint_vectors,
                              entry_scores::Vector{<:Real}, vis_threshold::Int)
    v_score = Dict{Int, Float64}()
    for endpoints in endpoint_vectors, (v, s) in zip(endpoints, entry_scores)
        v_score[v] = max(get(v_score, v, -Inf), s)
    end
    V_sorted = sort(V, by=v -> v_score[v], rev=true)
    return Set(V_sorted[1:min(vis_threshold, end)])
end

"""
    select_variable_edges(entry_scores, I, J, vis_threshold; bipartite=false, scored_on=:rows)
        -> (I_plot, J_plot, selected_indices)

Determine which edges of a graph, given by the vertex pairs `(I[k], J[k])`, to display, with
optional thresholding.

If `bipartite=false`, `I` and `J` are the source and destination vertices of the edges, both
subsets of one shared vertex set. If the number of vertices with at least one edge exceeds
`vis_threshold`, select an induced subgraph on the `vis_threshold` highest-scoring vertices,
where each vertex is scored by the maximum score of its incident edges (as either source or
destination). For an undirected graph, give each edge in only one orientation, `(i, j)` or
`(j, i)`; otherwise both are displayed.

If `bipartite=true`, `I` and `J` are subsets of the two separate partitions of a bipartite graph
(e.g. the rows and columns of a rectangular matrix), so `I[k]` and `J[k]` index different vertex
sets even when equal. Only the rows/columns given by `scored_on` (`:rows` for `I`, `:cols` for `J`)
is thresholded: if it has more than `vis_threshold` vertices with at least one edge, the
`vis_threshold` highest-scoring ones are kept, together with all of their edges. The other
dimension is not limited.

Returns the coordinates `I_plot`, `J_plot` of the entries to display, together with
`selected_indices`, their positions in the original `I`/`J` (so `I_plot == I[selected_indices]`
and each entry's data is row `selected_indices[k]` of every solver's `var_data`). Both are
needed downstream: `(I_plot[k], J_plot[k])` places panel `k` on the grid, while
`selected_indices[k]` locates its data.

If the entries selected for display contain repeated exact `(i, j)` coordinates (which can
happen when the variable encodes a multi-edge graph), only the highest-scoring occurrence of
each coordinate is kept and a warning is emitted per colliding coordinate (see
`dedup_edges_by_scores`).
"""
function select_variable_edges(entry_scores::Vector{<:Real},
                                I::Vector{Int}, J::Vector{Int},
                                vis_threshold::Int; bipartite::Bool=false,
                                scored_on::Symbol=:rows)
    vis_threshold > 0 || throw(ArgumentError("vis_threshold must be positive"))
    all(I .>= 1) || throw(ArgumentError("All I indices must be >= 1"))
    all(J .>= 1) || throw(ArgumentError("All J indices must be >= 1"))
    scored_on in (:rows, :cols) ||
        throw(ArgumentError("scored_on must be :rows or :cols, got :$scored_on"))

    if bipartite
        # Only the vertices on the scored partition/dimension are thresholded.
        scored_endpoints = scored_on == :rows ? I : J
        endpoint_vectors = (scored_endpoints,)
        is_selected = (subgraph, k) -> scored_endpoints[k] ∈ subgraph
    else
        # Each vertex is scored over the edges incident to it as either endpoint.
        endpoint_vectors = (I, J)
        is_selected = (subgraph, k) -> I[k] ∈ subgraph && J[k] ∈ subgraph
    end
    # Set of thresholded vertices with at least one edge
    V = union(unique.(endpoint_vectors)...)

    if length(V) <= vis_threshold
        # No need to filter. Collapse any repeated edges.
        I_plot, J_plot, selected_indices =
            dedup_edges_by_scores(I, J, collect(1:length(I)), entry_scores)
        return I_plot, J_plot, selected_indices
    end

    subgraph = top_scoring_vertices(V, endpoint_vectors, entry_scores, vis_threshold)
    # Keep the edges whose thresholded endpoints are in the set of selected vertices.
    # These indices point into the original I/J.
    selected_indices = findall(k -> is_selected(subgraph, k), 1:length(I))
    isempty(selected_indices) && throw(ArgumentError("No entries selected for display; consider increasing vis_threshold"))
    # Collapse any duplicate (i, j) coordinates among the selected entries.
    I_plot, J_plot, selected_indices = dedup_edges_by_scores(
        I[selected_indices], J[selected_indices],
        selected_indices, entry_scores[selected_indices]
    )
    return I_plot, J_plot, selected_indices
end
