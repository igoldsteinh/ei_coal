# this file is for creating coalescence trees from compartmental models
# with known infection histories
# and assuming that the relevant compartments are E and I (EI or SEIR both work)
using DataFrames
using Random
""" 
find_pair_coal_times(pair_ids, individ_frame)
find the coalescent time for a pair of individuals
# Arguments
-pair_ids: the pair of individual ids to find the coalescent time
-individ_frame: a dataframe containing infection, infectious, recovered times and transmission histories

"""
function find_pair_coal_times(pair_ids, individ_frame)
    # find the most recent common ancestor of the pair
    pair1_history = vcat(individ_frame.history[pair_ids[1]], pair_ids[1])
    pair2_history = vcat(individ_frame.history[pair_ids[2]], pair_ids[2])
    mrca = last(intersect(pair1_history, pair2_history))
    # two cases, first case, the mrca is one of the pair (i.e. one of the two individuals infected the other)
    # in this case, the coalescence time is the time of infection of the sampled individual who is not the mrca
    if mrca == pair_ids[1]
        mrca_index2 = findfirst(pair2_history .== mrca)
        coal_time = individ_frame.infec_time[pair2_history[mrca_index2 + 1]]
    elseif mrca == pair_ids[2]
        mrca_index1 = findfirst(pair1_history .== mrca)
        coal_time = individ_frame.infec_time[pair1_history[mrca_index1 + 1]]
    else 
        # now the third case is more complicated
        # in the third case the mrca is not one of the pair members
        # in which case, the mrca is some kind of ancestor of both pair members, 
        # for a single pair member, the mrca either infected them, or one of their more direct ancestors
        # in either case, the coalescence time is the minimum of the infection time of the lineage directly after the mrca
        # for each pair member 
        mrca_index1 = findfirst(pair1_history .== mrca)
        mrca_index2 = findfirst(pair2_history .== mrca)
        inf_time1 = individ_frame.infec_time[pair1_history[mrca_index1 + 1]]
        inf_time2 = individ_frame.infec_time[pair2_history[mrca_index2 + 1]]
        coal_time = minimum([inf_time1, inf_time2])
    end 
    return coal_time, mrca
end
"""
speed up find_pair_coal_times (modestly)
"""
function find_pair_coal_times2(pair_ids, individ_frame)
    # find the most recent common ancestor of the pair
    pair1_history = vcat(individ_frame.history[pair_ids[1]], pair_ids[1])
    pair2_history = vcat(individ_frame.history[pair_ids[2]], pair_ids[2])
    mrca = last(intersect(pair1_history, pair2_history))
    # two cases, first case, the mrca is one of the pair (i.e. one of the two individuals infected the other)
    # in this case, the coalescence time is the time of infection of the sampled individual who is not the mrca
    if mrca == pair_ids[1]
        mrca_index2 = findfirst(pair2_history .== mrca)
        coal_time = individ_frame.infec_time[pair2_history[mrca_index2 + 1]]
    elseif mrca == pair_ids[2]
        mrca_index1 = findfirst(pair1_history .== mrca)
        coal_time = individ_frame.infec_time[pair1_history[mrca_index1 + 1]]
    else 
        # now the third case is more complicated
        # in the third case the mrca is not one of the pair members
        # in which case, the mrca is some kind of ancestor of both pair members, 
        # for a single pair member, the mrca either infected them, or one of their more direct ancestors
        # in either case, the coalescence time is the minimum of the infection time of the lineage directly after the mrca
        # for each pair member 
        mrca_index1 = findfirst(pair1_history .== mrca)
        mrca_index2 = findfirst(pair2_history .== mrca)
        inf_time1 = individ_frame.infec_time[pair1_history[mrca_index1 + 1]]
        inf_time2 = individ_frame.infec_time[pair2_history[mrca_index2 + 1]]
        coal_time = minimum([inf_time1, inf_time2])
    end 
    return coal_time, mrca
end

""" 
find_samp_coal_times(samp_ids, individ_frame)
find the coalescence times for a set of sampled individuals
# Arguments
-samp_ids: list of individual ids to find all coalescent times of
-individ_frame: data frame containing infection,infectious,recovered times and infection histories

"""
function find_samp_coal_times(samp_ids, individ_frame)
    # for all pairs of samp_ids, find the coalescence time
    active_ids = samp_ids
    # create vector of length n-1 to store the coalescence times
    coalescence_times = Vector{Union{Missing, Float64}}(undef, length(samp_ids) -1)
    coal_pairs = []
    mrca_id = Vector{Int64}(undef, length(samp_ids) -1)
    for i in 1:(length(samp_ids) -1) 
    my_pairs = collect(combinations(active_ids, 2))
    pair_coal_times = Array{Union{Missing, Float64}}(undef, length(my_pairs), 2)
    for (j, pair) in enumerate(my_pairs)
        # print(pair)
        res = find_pair_coal_times(pair, individ_frame)
        pair_coal_times[j,1] = res[1]
        pair_coal_times[j,2] = res[2]
    end
    # print("active_ids")
    # print(active_ids)
    # find the maximum coalescence time (whichever one is sooner happens first)
    max_coal_time = maximum(pair_coal_times[:,1])
    # find the pair with the maximum coalescence time
    max_coal_pair = my_pairs[argmax(pair_coal_times[:,1])]
    # find the mrca of the pair with the max coalescence time
    max_coal_mrca = pair_coal_times[argmax(pair_coal_times[:,1]),2]
    # remove the pair with the maximum coalescence time from the active_ids
    active_ids = setdiff(active_ids, max_coal_pair)
    old_active_ids = copy(active_ids)
    # now we need to find the ancestors of all of the remaining lineages at the time of the coalescence
    for old_lineage in old_active_ids
        # find the vector of possible lineages, which includes the active_id
        possible_lineages = vcat(individ_frame.history[old_lineage], old_lineage)
        # filter individ_frame for the possible lineages
        possible_frame = filter(row -> row.id in possible_lineages, eachrow(individ_frame))
        possible_frame = filter(row -> row.infec_time < max_coal_time, possible_frame)
        possible_frame = filter(row -> ismissing(row.recov_time) || row.recov_time > max_coal_time, possible_frame)
        # find all lineages who were infected before the coalescence time and who have not yet recovered by the coalescence time
        # eligible_lineages = filter(row -> row.infec_time < max_coal_time && (ismissing(row.recov_time) || row.recov_time > max_coal_time) && row.id in possible_lineages, eachrow(individ_frame))
        # pick the most recent lineage 
        # print("old lineage")
        # print(old_lineage)
        closest_lineage = possible_frame.id[argmax(possible_frame.infec_time)]
        # remove the active_id from the active_ids
        active_ids = filter!(x -> x != old_lineage, active_ids)
        # add the closest lineage to the active_ids
        push!(active_ids, closest_lineage)
    end

    # add the mrca to the active_ids
    push!(active_ids, max_coal_mrca)
    # store the minimum coalescence time
    coalescence_times[i] = max_coal_time
    push!(coal_pairs, max_coal_pair)
    mrca_id[i] = max_coal_mrca
    # print(i)
    # print("active_ids")
    # print(active_ids)
    end
    # return the coalescence times
    return coalescence_times, coal_pairs, mrca_id
end
"""
rewrite find_samp_coal_times to make it faster
"""
function find_samp_coal_times2(samp_ids, individ_frame)
    # for all pairs of samp_ids, find the coalescence time
    active_ids = samp_ids
    # create vector of length n-1 to store the coalescence times
    coalescence_times = Vector{Union{Missing, Float64}}(undef, length(samp_ids) -1)
    coal_pairs = []
    mrca_id = Vector{Int64}(undef, length(samp_ids) -1)
    for i in 1:(length(samp_ids) -1) 
    my_pairs = collect(combinations(active_ids, 2))
    pair_coal_times = Array{Union{Missing, Float64}}(undef, length(my_pairs), 2)
    for (j, pair) in enumerate(my_pairs)
        # print(pair)
        res = find_pair_coal_times(pair, individ_frame)
        pair_coal_times[j,1] = res[1]
        pair_coal_times[j,2] = res[2]
    end
    # print("active_ids")
    # print(active_ids)
    # find the maximum coalescence time (whichever one is sooner happens first)
    max_coal_time = maximum(pair_coal_times[:,1])
    # find the pair with the maximum coalescence time
    max_coal_pair = my_pairs[argmax(pair_coal_times[:,1])]
    # find the mrca of the pair with the max coalescence time
    max_coal_mrca = pair_coal_times[argmax(pair_coal_times[:,1]),2]
    # remove the pair with the maximum coalescence time from the active_ids
    active_ids = setdiff(active_ids, max_coal_pair)
    old_active_ids = copy(active_ids)
    # now we need to find the ancestors of all of the remaining lineages at the time of the coalescence
    for old_lineage in old_active_ids
        # find the vector of possible lineages, which includes the active_id
        possible_lineages = vcat(individ_frame.history[old_lineage], old_lineage)
        # filter individ_frame for the possible lineages
        possible_frame = filter(row -> row.id in possible_lineages &&
         row.infec_time < max_coal_time && 
         (ismissing(row.recov_time) || row.recov_time > max_coal_time), eachrow(individ_frame))
        # find all lineages who were infected before the coalescence time and who have not yet recovered by the coalescence time
        # eligible_lineages = filter(row -> row.infec_time < max_coal_time && (ismissing(row.recov_time) || row.recov_time > max_coal_time) && row.id in possible_lineages, eachrow(individ_frame))
        # pick the most recent lineage 
        # print("old lineage")
        # print(old_lineage)
        closest_lineage = possible_frame.id[argmax(possible_frame.infec_time)]
        # remove the active_id from the active_ids
        active_ids = filter!(x -> x != old_lineage, active_ids)
        # add the closest lineage to the active_ids
        push!(active_ids, closest_lineage)
    end

    # add the mrca to the active_ids
    push!(active_ids, max_coal_mrca)
    # store the minimum coalescence time
    coalescence_times[i] = max_coal_time
    push!(coal_pairs, max_coal_pair)
    mrca_id[i] = max_coal_mrca
    # print(i)
    # print("active_ids")
    # print(active_ids)
    end
    # return the coalescence times
    return coalescence_times, coal_pairs, mrca_id
end
"""
find_soonest_coal_times(samp_ids, individ_frame)
find the soonest coalescence times for a set of sampled individuals
this should be faster than finding all of them 
"""
function find_soonest_coal_times(samp_ids, individ_frame)
    # for all pairs of samp_ids, find the coalescence time
    active_ids = samp_ids
    # create vector of length n-1 to store the coalescence times
    my_pairs = collect(combinations(active_ids, 2))
    pair_coal_times = Array{Union{Missing, Float64}}(undef, length(my_pairs), 2)
    for (j, pair) in enumerate(my_pairs)
        # print(pair)
        res = find_pair_coal_times(pair, individ_frame)
        pair_coal_times[j,1] = res[1]
        pair_coal_times[j,2] = res[2]
    end
    max_coal_time = maximum(pair_coal_times[:,1])
    # find the pair with the maximum coalescence time
    max_coal_pair = my_pairs[argmax(pair_coal_times[:,1])]
    # find the mrca of the pair with the max coalescence time
    max_coal_mrca = pair_coal_times[argmax(pair_coal_times[:,1]),2]
    return max_coal_time, max_coal_pair, max_coal_mrca
end
""" 
find_soonest_i2e_time(I_ids, individ_frame)
find the first time when a lineage switches from I to E 
# Arguments
I_ids: list of ids in the I state 
individ_frame: frame containing infectious times and transmission histories
"""
function find_soonest_i2e_time(I_ids, individ_frame)
    # find the soonest time when the lineages switch from I to E 
    #this is simply the time at which they became infectious
    infectious_times = individ_frame.infectious_time[I_ids]
    # find the maximum time
    soonest_time = maximum(infectious_times)
    # find the id of the individual who became infectious at the soonest time
    soonest_id = I_ids[argmax(infectious_times)]
    return soonest_time, soonest_id
end 
""" 
find_soonest_e2i_time(E_ids, individ_frame)
find the first time when a lineage switches from E to I
# Arguments
E_ids: list of ids in the E state 
individ_frame: data frame containing infection times and transmission histories
"""
function find_soonest_e2i_time(E_ids, individ_frame)
    # find the soonest time when the lineages switch from E to I
    # this is the time at which they were infected
    inf_times = individ_frame.infec_time[E_ids]
    # find the minimum time
    soonest_time = maximum(inf_times)
    # find the id of the individual who was infected at the soonest time
    soonest_id = E_ids[argmax(inf_times)]
    # find their infector
    infector = individ_frame.history[soonest_id][end]
    return soonest_time, soonest_id, infector
end
""" 
create_ei_coal_tree(samp_ids, individ_frame, samp_time)
create a coalescent tree from a set of sampled individuals
output is a dataframe with columns with times, event types and counts of lineages
this will be done in forward time for sanity's sake
reverse time will be calculated separately
# Arguments
samp_ids: list of ids to construct the tree from
individ_frame: data frame containing all infection,infectious,recovered times and transmission histories
samp_time: the time at which sampling occurred (in forward time)
"""
function create_ei_coal_tree(samp_ids, individ_frame, samp_time)
    # create a dataframe to store the tree
    # samp_ids = [342,351,53,251,345,207,214,170,292,270]
    # find all of samp_ids which are descendants of 53
    tree = DataFrame(time = Float64[], event = String[], nE = Int64[], nI = Int64[], samp_ids = Any[], og_ids = Any[])
    # add time 0 in 
    t = samp_time
    og_ids = DataFrame(samp_ids = samp_ids, og_ids = samp_ids)

    push!(tree, (t, "start", 0, length(samp_ids), samp_ids, og_ids))
    I_ids = samp_ids
    E_ids = []
    tree_coal_pairs = []
    # find the first time when a lineage switches from I to E
    while length(samp_ids) > 1
        # calculate the soonest i2e, e2i and coalescence times
        if length(I_ids) > 0
            i2e_time, i2e_id = find_soonest_i2e_time(I_ids, individ_frame)
        else
            i2e_time = -1
        end
        if length(E_ids) > 0
            e2i_time, e2i_id, infector_id = find_soonest_e2i_time(E_ids, individ_frame)
        else
            e2i_time = -1
        end
        # find the coalescence times
        coalescence_times, coal_pairs, mrca_id = find_samp_coal_times(samp_ids, individ_frame)
        max_coal_time = maximum(coalescence_times)
        max_coal_pair = coal_pairs[argmax(coalescence_times)]
        max_mrca_id = mrca_id[argmax(coalescence_times)]
        # execute the soonest event
        next_time = maximum([i2e_time, e2i_time, max_coal_time])
        if next_time == i2e_time
            # remove the id from the I_ids and add it to the E_ids
            push!(E_ids, i2e_id)
            I_ids = filter(x -> x != i2e_id, I_ids)
            # update time to the i2e_time
            t = i2e_time
            # find the og_id from the i2e_id 
            og_id = og_ids[og_ids.samp_ids .== i2e_id, :og_ids][1]

            push!(tree, (t, "I2E", length(E_ids), length(I_ids), samp_ids, og_id))
        elseif next_time == e2i_time && next_time != max_coal_time
            # find the og_id from the e2i_id
            og_id = og_ids[og_ids.samp_ids .== e2i_id, :og_ids][1]
            # update the og_ids to include the infector id
            og_ids[og_ids.samp_ids .== e2i_id, :samp_ids] = infector_id
            # remove the id from the E_ids also from samp_ids
            E_ids = filter(x -> x != e2i_id, E_ids)
            samp_ids = filter(x -> x != e2i_id, samp_ids)
            # add the infector id to the I_ids and to samp_ids
            push!(I_ids, infector_id)
            push!(samp_ids, infector_id)
            # update time to the e2i_time
            t = e2i_time
            push!(tree, (t, "E2I", length(E_ids), length(I_ids), samp_ids, og_id))
        else
            # remove the pair with the minimum coalescence time from the samp_ids and from E and I ids
            # find the og_ids for the max_coal_pair
            og_ids_pair = og_ids[og_ids.samp_ids in max_coal_pair, :og_ids]
            # max_mrca_id needs to be added to the og_ids
            # its not entirely clear how that should be done
            # ah, its whoever was in state I at the time of coalescence 
            current_i_id = filter(x -> x ∈ max_coal_pair, I_ids)
            og_ids[og_ids.samp_ids .== current_i_id, :samp_ids] = max_mrca_id
            samp_ids = setdiff(samp_ids, max_coal_pair)
            E_ids = filter(x -> x ∉ max_coal_pair, E_ids)
            I_ids = filter(x -> x ∉ max_coal_pair, I_ids)
            # add the mrca to the I_ids and to samp_ids
            push!(I_ids, max_mrca_id)
            push!(samp_ids, max_mrca_id)
            push!(tree_coal_pairs, max_coal_pair)
            # update time to the min_coal_time
            t = max_coal_time
            push!(tree, (t, "coal", length(E_ids), length(I_ids), samp_ids, og_ids_pair))
        end
    end
    return tree, tree_coal_pairs
end

""" 
create_ei_coal_tree_hetchron(samp_ids, individ_frame, samp_time)
create a coalescent tree from a set of sampled individuals
output is a dataframe with columns with times, event types and counts of lineages
this will be done in forward time for sanity's sake
reverse time will be calculated separately
# Arguments
samp_ids: list of ids to construct the tree from
individ_frame: data frame containing all infection,infectious,recovered times and transmission histories
samp_time: the time at which sampling occurred (in forward time)
forward_samp_times: the times at which the forward samples were taken
forward_samp_ids: the ids of the individuals sampled at the forward sample times
the rule is that after you've been sampled, none of your descendendants created post sampling can be sampled
"""
function create_ei_coal_tree_hetchron(last_samp_ids, last_samp_time, forward_samp_times, forward_samp_ids, individ_frame)
    # create a dataframe to store the tree
    # samp_ids = [342,351,53,251,345,207,214,170,292,270]
    # find all of samp_ids which are descendants of 53
    time = Float64[]
    event = String[]
    nE = Int64[]
    nI = Int64[]
    samp_ids_overtime = Vector{Vector{Int64}}()
    og_ids_overtime = Vector{Vector{Int64}}()
    # add time 0 in 
    t = last_samp_time
    samp_ids = zeros(Int64, length(last_samp_ids))
    samp_ids .= last_samp_ids
    active_samp_ids = zeros(Int64, length(last_samp_ids))
    active_samp_ids .= last_samp_ids
    og_ids = zeros(Int64, length(last_samp_ids))
    og_ids .= last_samp_ids
    next_samp_idx = length(forward_samp_times)
    next_samp_time = forward_samp_times[next_samp_idx]
    next_samp_ids = forward_samp_ids[next_samp_idx]
    push!(time, t)
    push!(event, "start")
    push!(nE, 0)
    push!(nI, length(samp_ids))
    push!(samp_ids_overtime, deepcopy(active_samp_ids))
    push!(og_ids_overtime, og_ids)
    I_ids = samp_ids
    E_ids = []
    tree_coal_pairs = []
    # find the first time when a lineage switches from I to E
    while length(active_samp_ids) > 1
        # calculate the soonest i2e, e2i and coalescence times
        if length(I_ids) > 0
            i2e_time, i2e_id = find_soonest_i2e_time(I_ids, individ_frame)
        else
            i2e_time = -1
        end
        if length(E_ids) > 0
            e2i_time, e2i_id, infector_id = find_soonest_e2i_time(E_ids, individ_frame)
        else
            e2i_time = -1
        end
        # find the coalescence times
        coalescence_times, coal_pairs, mrca_id = find_samp_coal_times(active_samp_ids, individ_frame)
        max_coal_time = maximum(coalescence_times)
        max_coal_pair = coal_pairs[argmax(coalescence_times)]
        max_mrca_id = mrca_id[argmax(coalescence_times)]
        # execute the soonest event
        next_time = maximum([i2e_time, e2i_time, max_coal_time, next_samp_time])
        if next_time == next_samp_time 
            # sample new lineages from the available pool 
            samp_ids = vcat(samp_ids, next_samp_ids)
            active_samp_ids = vcat(active_samp_ids, next_samp_ids)
            og_ids = vcat(og_ids, next_samp_ids)
            I_ids = vcat(I_ids, next_samp_ids)
            # update time to the next_samp_time
            t = next_samp_time
            push!(time, t)
            push!(event, "samp")
            push!(nE, length(E_ids))
            push!(nI, length(I_ids))
            push!(samp_ids_overtime, deepcopy(active_samp_ids))
            push!(og_ids_overtime, next_samp_ids)

            # update the next_samp_idx
            if next_samp_idx > 1
                next_samp_idx = next_samp_idx - 1
                next_samp_time = forward_samp_times[next_samp_idx]
                next_samp_ids = forward_samp_ids[next_samp_idx]
    
            else
                next_samp_idx = next_samp_idx
                next_samp_time = -1
                next_samp_ids = []
            end
        elseif next_time == i2e_time && next_time != max_coal_time
            # remove the id from the I_ids and add it to the E_ids
            push!(E_ids, i2e_id)
            I_ids = filter(x -> x != i2e_id, I_ids)
            # update time to the i2e_time
            t = i2e_time
            push!(time, t)
            push!(event, "I2E")
            push!(nE, length(E_ids))
            push!(nI, length(I_ids))
            push!(samp_ids_overtime, deepcopy(active_samp_ids))
            # find the og_id from the i2e_id 
            i2e_index = findfirst(samp_ids .== i2e_id)
            current_og_id = og_ids[i2e_index]
            push!(og_ids_overtime, [current_og_id])

        elseif next_time == e2i_time && next_time != max_coal_time
            # find the og_id from the e2i_id
            e2i_index = findfirst(samp_ids .== e2i_id)
            current_og_id = og_ids[e2i_index]
            # remove the id from the E_ids also from samp_ids
            E_ids = filter(x -> x != e2i_id, E_ids)
            samp_ids[findfirst(samp_ids .== e2i_id)] = infector_id
            active_samp_ids[findfirst(active_samp_ids .== e2i_id)] = infector_id
            # add the infector id to the I_ids and to samp_ids
            push!(I_ids, infector_id)
            # update time to the e2i_time
            t = e2i_time
            push!(time, t)
            push!(event, "E2I")
            push!(nE, length(E_ids))
            push!(nI, length(I_ids))
            push!(samp_ids_overtime, deepcopy(active_samp_ids))
            push!(og_ids_overtime, [current_og_id])
        else
            # update og_ids and samp_ids 
            current_og_ids = og_ids[findall(x -> x in max_coal_pair, samp_ids)]
            current_i_id = filter(x -> x ∈ max_coal_pair, I_ids)
            i_id_index = findfirst(samp_ids .== current_i_id)  
            samp_ids[i_id_index] = max_mrca_id

            # remove the pair with the minimum coalescence time from the samp_ids and from E and I ids
            active_samp_ids = setdiff(active_samp_ids, max_coal_pair)
            E_ids = filter(x -> x ∉ max_coal_pair, E_ids)
            I_ids = filter(x -> x ∉ max_coal_pair, I_ids)
            # add the mrca to the I_ids and to samp_ids
            push!(I_ids, max_mrca_id)
            push!(active_samp_ids, max_mrca_id)
            push!(tree_coal_pairs, max_coal_pair)
            # update time to the min_coal_time
            t = max_coal_time
            push!(time, t)
            push!(event, "coal")
            push!(nE, length(E_ids))
            push!(nI, length(I_ids))
            push!(samp_ids_overtime, deepcopy(active_samp_ids))
            push!(og_ids_overtime, current_og_ids)
        end
    end
    return DataFrame(time=time, event=event, nE=nE, nI=nI, samp_ids=samp_ids_overtime, og_ids = og_ids_overtime), tree_coal_pairs
end
""" 
create_ei_coal_tree_hetchron2(samp_ids, individ_frame, samp_time)
trying to make this faster, don't want to break what we already have
create a coalescent tree from a set of sampled individuals
output is a dataframe with columns with times, event types and counts of lineages
this will be done in forward time for sanity's sake
reverse time will be calculated separately
# Arguments
samp_ids: list of ids to construct the tree from
individ_frame: data frame containing all infection,infectious,recovered times and transmission histories
samp_time: the time at which sampling occurred (in forward time)
forward_samp_times: the times at which the forward samples were taken
forward_samp_ids: the ids of the individuals sampled at the forward sample times
the rule is that after you've been sampled, none of your descendendants created post sampling can be sampled
"""
function create_ei_coal_tree_hetchron2(last_samp_ids, last_samp_time, forward_samp_times, forward_samp_ids, individ_frame)
    # create a dataframe to store the tree
    # samp_ids = [342,351,53,251,345,207,214,170,292,270]
    # find all of samp_ids which are descendants of 53
    time = Float64[]
    event = String[]
    nE = Int64[]
    nI = Int64[]
    samp_ids_overtime = Vector{Vector{Int64}}()
    og_ids_overtime = Vector{Vector{Int64}}()
    # add time 0 in 
    t = last_samp_time
    samp_ids = zeros(Int64, length(last_samp_ids))
    samp_ids .= last_samp_ids
    active_samp_ids = zeros(Int64, length(last_samp_ids))
    active_samp_ids .= last_samp_ids
    og_ids = zeros(Int64, length(last_samp_ids))
    og_ids .= last_samp_ids
    next_samp_idx = length(forward_samp_times)
    if length(forward_samp_times) == 0
        next_samp_time = -1
        next_samp_ids = []
    else
        next_samp_time = forward_samp_times[next_samp_idx]
        next_samp_ids = forward_samp_ids[next_samp_idx]
    end
    push!(time, t)
    push!(event, "samp")
    push!(nE, 0)
    push!(nI, length(samp_ids))
    push!(samp_ids_overtime, deepcopy(active_samp_ids))
    push!(og_ids_overtime, og_ids)
    I_ids = samp_ids
    E_ids = []
    tree_coal_pairs = []
    max_coal_time = -1
    max_coal_pair = []
    max_mrca_id = -1
    # find the first time when a lineage switches from I to E
    while length(active_samp_ids) > 1
        # calculate the soonest i2e, e2i and coalescence times
        if length(I_ids) > 0
            i2e_time, i2e_id = find_soonest_i2e_time(I_ids, individ_frame)
        else
            i2e_time = -1
        end
        if length(E_ids) > 0
            e2i_time, e2i_id, infector_id = find_soonest_e2i_time(E_ids, individ_frame)
        else
            e2i_time = -1
        end
        # find the coalescence times
        # this only needs to be done if the previous event was a sampling event
        # or if it was a coalescence event
        # or if its an e2i event
        if event[end] == "samp" || event[end] == "coal" || event[end] == "E2I"
            # find the coalescence times
            max_coal_time, max_coal_pair, max_mrca_id = find_soonest_coal_times(active_samp_ids, individ_frame)
        end
        # execute the soonest event
        next_time = maximum([i2e_time, e2i_time, max_coal_time, next_samp_time])
        if next_time == next_samp_time 
            # sample new lineages from the available pool 
            samp_ids = vcat(samp_ids, next_samp_ids)
            active_samp_ids = vcat(active_samp_ids, next_samp_ids)
            og_ids = vcat(og_ids, next_samp_ids)
            I_ids = vcat(I_ids, next_samp_ids)
            # update time to the next_samp_time
            t = next_samp_time
            push!(time, t)
            push!(event, "samp")
            push!(nE, length(E_ids))
            push!(nI, length(I_ids))
            push!(samp_ids_overtime, deepcopy(active_samp_ids))
            push!(og_ids_overtime, next_samp_ids)

            # update the next_samp_idx
            if next_samp_idx > 1
                next_samp_idx = next_samp_idx - 1
                next_samp_time = forward_samp_times[next_samp_idx]
                next_samp_ids = forward_samp_ids[next_samp_idx]
    
            else
                next_samp_idx = next_samp_idx
                next_samp_time = -1
                next_samp_ids = []
            end
        elseif next_time == i2e_time && next_time != max_coal_time
            # remove the id from the I_ids and add it to the E_ids
            push!(E_ids, i2e_id)
            I_ids = filter(x -> x != i2e_id, I_ids)
            # update time to the i2e_time
            t = i2e_time
            push!(time, t)
            push!(event, "I2E")
            push!(nE, length(E_ids))
            push!(nI, length(I_ids))
            push!(samp_ids_overtime, deepcopy(active_samp_ids))
            # find the og_id from the i2e_id 
            i2e_index = findfirst(samp_ids .== i2e_id)
            current_og_id = og_ids[i2e_index]
            push!(og_ids_overtime, [current_og_id])

        elseif next_time == e2i_time && next_time != max_coal_time
            # find the og_id from the e2i_id
            e2i_index = findfirst(samp_ids .== e2i_id)
            current_og_id = og_ids[e2i_index]
            # remove the id from the E_ids also from samp_ids
            E_ids = filter(x -> x != e2i_id, E_ids)
            samp_ids[findfirst(samp_ids .== e2i_id)] = infector_id
            active_samp_ids[findfirst(active_samp_ids .== e2i_id)] = infector_id
            # add the infector id to the I_ids and to samp_ids
            push!(I_ids, infector_id)
            # update time to the e2i_time
            t = e2i_time
            push!(time, t)
            push!(event, "E2I")
            push!(nE, length(E_ids))
            push!(nI, length(I_ids))
            push!(samp_ids_overtime, deepcopy(active_samp_ids))
            push!(og_ids_overtime, [current_og_id])
        else
            # update og_ids and samp_ids 
            current_og_ids = og_ids[findall(x -> x in max_coal_pair, samp_ids)]
            current_i_id = filter(x -> x ∈ max_coal_pair, I_ids)
            i_id_index = findfirst(samp_ids .== current_i_id)  
            samp_ids[i_id_index] = max_mrca_id

            # remove the pair with the minimum coalescence time from the samp_ids and from E and I ids
            active_samp_ids = setdiff(active_samp_ids, max_coal_pair)
            E_ids = filter(x -> x ∉ max_coal_pair, E_ids)
            I_ids = filter(x -> x ∉ max_coal_pair, I_ids)
            # add the mrca to the I_ids and to samp_ids
            push!(I_ids, max_mrca_id)
            push!(active_samp_ids, max_mrca_id)
            push!(tree_coal_pairs, max_coal_pair)
            # update time to the min_coal_time
            t = max_coal_time
            push!(time, t)
            push!(event, "coal")
            push!(nE, length(E_ids))
            push!(nI, length(I_ids))
            push!(samp_ids_overtime, deepcopy(active_samp_ids))
            push!(og_ids_overtime, current_og_ids)
        end
    end
    return DataFrame(time=time, event=event, nE=nE, nI=nI, samp_ids=samp_ids_overtime, og_ids = og_ids_overtime), tree_coal_pairs
end

"""
sample_infected_ids(individ_frame, last_samp_time, forward_samp_times, forward_samp_lin, num_samps, iso)
sample ids from an outbreak
to be used as input in tree creation
"""
function sample_infected_ids(individ_frame, last_samp_time, forward_samp_times, forward_samp_lin, num_samps, iso)
    # first check if isochronous or not 
    if iso == true 
        # sample anyone who is currently infectious at the last_samp_time
        individ_frame_clean = dropmissing(individ_frame, [:infectious_time])
        eligible_individuals = filter(row -> row.infectious_time < last_samp_time && (ismissing(row.recov_time) || row.recov_time > last_samp_time), eachrow(individ_frame_clean))
        samp_ids = [sample(eligible_individuals.id, num_samps, replace=false)]
    else 
        # iso is false, that is we are doing heterechronous sampling 
        forward_samp_times = sort(forward_samp_times)
        forward_samp_ids_flat_vec = []
        samp_ids = []
        prev_samp_time = 0
        individ_frame_clean = dropmissing(individ_frame, [:infectious_time])
        for t in 1:length(forward_samp_times) 
            samp_time = forward_samp_times[t]
            samp_lin = forward_samp_lin[t]
            # sample from the individuals who are infectious at the time of sampling
            # and who are not recovered at the time of sampling
            # and who cannot be infected by a previously sampled individual after the time of sampling 
            # Filter the DataFrame directly using column access
            individ_frame_clean = filter(row -> !(row[:infec_time] > prev_samp_time && any(id -> id in forward_samp_ids_flat_vec, row[:history])) && 
            !(row[:id] in forward_samp_ids_flat_vec), individ_frame_clean)
            eligible_individuals = filter(row -> row[:infectious_time] < samp_time && 
                (ismissing(row[:recov_time]) || row[:recov_time] > samp_time), individ_frame_clean)

            new_samp_ids = sample(eligible_individuals.id, samp_lin, replace=false)
            push!(samp_ids, new_samp_ids)
            forward_samp_ids_flat_vec = vcat(forward_samp_ids_flat_vec, new_samp_ids)
            prev_samp_time = samp_time
            test_frame = filter(row -> new_samp_ids[1] in row[:history], individ_frame_clean)
            testing = filter(row -> 219 in row[:history], individ_frame_clean)
        end
        # filter for individuals who are descended from forward_samp_ids_flat_vec
        individ_frame_clean = filter(row -> !(row[:infec_time] > prev_samp_time && any(id -> id in forward_samp_ids_flat_vec, row[:history])) && 
            !(row[:id] in forward_samp_ids_flat_vec), individ_frame_clean)
        eligible_individuals = filter(row -> row[:infectious_time] < last_samp_time && 
                (ismissing(row[:recov_time]) || row[:recov_time] > last_samp_time), individ_frame_clean)
        # sample from the individuals who are still eligible at the last sampling time
        last_samp_ids = sample(eligible_individuals.id, num_samps, replace=false)
        push!(samp_ids, last_samp_ids)
    end 
    return samp_ids
end
