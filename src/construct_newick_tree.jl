using DataFrames
using Random
function construct_newick_tree(df::DataFrame, reverse_samp_lin::Vector{Int})
    # Sort by time (increasing)
    df_sorted = sort(df, :time)
    # Each active lineage is a tuple: (Set of IDs, Newick label, time)
    active = Vector{Tuple{Set{Int}, String, Float64}}()
    leaf_count = 0
    samp_idx = 1
    for row in eachrow(df_sorted)
        if row.event == "samp"
            for i in 1:reverse_samp_lin[samp_idx]
                leaf_count += 1
                current_id = row.og_ids[i]
                push!(active, (Set([current_id]), "L$current_id", row.time))
            end
            samp_idx += 1
        elseif row.event == "coal"
            # Find the two active lineages that contain the og_ids to coalesce
            idxs = [findfirst(x -> row.og_ids[i] in x[1], active) for i in 1:2]
            idx1, idx2 = sort(idxs)  # sort so we remove higher index first
            set1, l1, t1 = active[idx2]
            set2, l2, t2 = active[idx1]
            # Remove the two selected lineages
            deleteat!(active, idx2)
            deleteat!(active, idx1)
            # Branch lengths: time since last event for each lineage
            bl1 = row.time - t1
            bl2 = row.time - t2
            # Newick string for this coalescence
            new_label = "($l1:$(bl1),$l2:$(bl2))"
            # New set of IDs
            new_set = union(set1, set2)
            # Push new lineage with current time
            push!(active, (new_set, new_label, row.time))
        end
    end
    # At the end, there should be one lineage left
    if length(active) != 1
        error("Final tree does not have exactly one root!")
    end
    return active[1][2] * ";"
end

