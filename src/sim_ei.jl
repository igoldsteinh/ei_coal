# simulate EI model 
using DataFrames
using Random
function sim_ei(I_init::Int64, alpha::Float64, nu::Float64, gamma::Float64, stop_time::Float64)
    # Setting up initial states and frames
    t = 0.0
    num_exposed = 0
    num_infectious = I_init
    num_recovered = 0
    rt = alpha/nu 
    state_frame = [t num_exposed num_infectious num_recovered alpha rt]

    # initialize ids and individ_frame
    exposed_ids = Int[]
    infectious_ids = Int[]
    recovered_ids = Int[]
    # Initialize the last id
    current_id = I_init
    
    # Tracking history and times
    individ_frame = DataFrame(id=Int[], 
                              infec_time=Union{Missing, Float64}[], 
                              infectious_time=Union{Missing, Float64}[], 
                              recov_time=Union{Missing, Float64}[], 
                              history=Any[])  # List of who infected them

    # Initialize initial exposed individuals
    for i in 1:I_init
        push!(infectious_ids, i)
        push!(individ_frame, (i, -1, 0, missing, [i]))
    end
    # simulate till end point reached
    while (num_exposed > 0 || num_infectious > 0) && t < stop_time
        # Time to the next event
        next_event = rand(Exponential(1 / (alpha * num_infectious + gamma * num_exposed + nu * num_infectious)))

        # Update time
        t += next_event

        # Choose which event happens proportional to the rates
        which_event = sample(["infection", "infectious", "recovery"], Weights([alpha * num_infectious, gamma * num_exposed, nu * num_infectious]))

        # If the event is infection
        if which_event == "infection"
            # update current id
            current_id += 1
            # choose an infector from among the infectious individuals at random
            infector_id = sample(infectious_ids)
            # Update states
            num_exposed += 1
            push!(exposed_ids, current_id)
            push!(individ_frame, (current_id, t, missing, missing, unique([individ_frame.history[infector_id]..., infector_id])))
        # If the event is becoming infectious
        elseif which_event == "infectious"
            # choose an exposed individual at random
            exposed_id = sample(exposed_ids)
            # Update states
            num_exposed -= 1
            num_infectious += 1
            push!(infectious_ids, exposed_id)
            # remove exposed_id from exposed_ids
            exposed_ids = filter!(x -> x != exposed_id, exposed_ids)

            individ_frame.infectious_time[exposed_id] = t
        # If the event is recovery
        else
            # choose an infectious individual at random
            infectious_id = sample(infectious_ids)
            # Update states
            num_infectious -= 1
            num_recovered += 1
            push!(recovered_ids, infectious_id)
            # remove infectious_id from infectious_ids
            infectious_ids = filter!(x -> x != infectious_id, infectious_ids)
            individ_frame.recov_time[infectious_id] = t
        end
        # update state_frame
        rt = alpha/nu
        state_frame = vcat(state_frame, [t num_exposed num_infectious num_recovered alpha rt])
    end
    return (individ_frame, DataFrame(state_frame, [:time, :E, :I, :R, :alpha, :rt]))
end


# sample individuals but only infectious ones
function samp_individs_nolatent(individ_frame::DataFrame, t::Float64, n::Int)
    # Filter eligible individuals who were infectious at time T
        # Drop rows with missing infectious_time 
        individ_frame_clean = dropmissing(individ_frame, [:infectious_time])

    eligible_individuals = filter(row -> row.infectious_time < t && (ismissing(row.recov_time) || row.recov_time > t), eachrow(individ_frame_clean))


    # Sample n individuals from the eligible individuals
    samp_ids = sample(eligible_individuals.id, n, replace=false)  # `replace=false` ensures sampling without replacement

    return samp_ids
end

function find_true_ancestors_latent(t::Float64, samp_ids::Vector{Int}, individ_frame::DataFrame)
    ancestor_frame = DataFrame(samp_id=samp_ids, 
    ancestor_id=Vector{Union{Missing, Int64}}(undef, length(samp_ids)), 
    ancestor_state=Vector{Union{Missing, String}}(undef, length(samp_ids)),
    time=t)

    # Calculate ancestors at time t for sampled individuals
    for (i, samp_id) in enumerate(samp_ids)
        infec_time = individ_frame.infec_time[samp_id]
        infectious_time = individ_frame.infectious_time[samp_id]

        if infec_time < t 
            # If the time of interest is after infection, you are your own ancestor
            ancestor_frame.ancestor_id[i] = samp_id
            # the state of the ancestor is the state of the individual at time t
        elseif infec_time >= t
            # If your infection happens after time t
            ancestor_history = individ_frame.history[samp_id]
            possible_individuals = filter(row -> row.id in ancestor_history, eachrow(individ_frame))
            eligible_ancestors = filter(row -> row.infec_time < t, possible_individuals)

            if size(eligible_ancestors, 1) > 0
                closest_ancestor = argmax(eligible_ancestors.infec_time)
                ancestor_frame.ancestor_id[i] = eligible_ancestors.id[closest_ancestor]
            end
        end
        # record the state of the ancestor at time t
        if coalesce(individ_frame.infectious_time[ancestor_frame.ancestor_id[i]] > t, true)
            ancestor_frame.ancestor_state[i] = "E"
        else
            ancestor_frame.ancestor_state[i] = "I"
        end
    end

    return ancestor_frame
end

function ancestry_over_time_latent(start_time::Float64, grid_size::Float64, samp_ids::Vector{Int}, individ_frame::DataFrame)
    num_ancestors = length(samp_ids)
    num_ancestors_frame = DataFrame(time=collect(0:grid_size:start_time), 
    num_ancestors=Vector{Union{Missing, Int64}}(undef,length(0:grid_size:start_time)),
    num_ancestors_e = Vector{Union{Missing, Int64}}(undef,length(0:grid_size:start_time)),
    num_ancestors_i = Vector{Union{Missing, Int64}}(undef,length(0:grid_size:start_time)))

    num_ancestors_frame.num_ancestors[num_ancestors_frame.time .== start_time] .= num_ancestors
    # count the number of e ancestors at the start time which is the number of samp_ids whose infectious time is after the start time or whose infectious time is missing
    num_ancestors_e = length(filter(id -> id in individ_frame.id[coalesce.(individ_frame.infectious_time .> start_time, true)], samp_ids))
    # find the ids for whom 
    num_ancestors_frame.num_ancestors_e[num_ancestors_frame.time .== start_time] .= num_ancestors_e
    # count the number of i ancestors at the start time which is the number of samp_ids whose infectious time is before the start time
    # this is simply the difference between the total ancestors and the number of i1 ancestors
    num_ancestors_i = num_ancestors - num_ancestors_e
    num_ancestors_frame.num_ancestors_i[num_ancestors_frame.time .== start_time] .= num_ancestors_i

    current_time = start_time

    while num_ancestors > 1 && current_time > 0
        current_time -= grid_size
        ancestors = find_true_ancestors_latent(current_time, samp_ids, individ_frame)
        num_ancestors = length(unique(ancestors.ancestor_id))
        #count the number of i1 ancestors which are unique
        num_ancestors_e = length(unique(ancestors.ancestor_id[ancestors.ancestor_state .== "E"]))
        #count the number of i2 ancestors which are unique
        num_ancestors_i = length(unique(ancestors.ancestor_id[ancestors.ancestor_state .== "I"]))
        num_ancestors_frame.num_ancestors[num_ancestors_frame.time .== current_time] .= num_ancestors
        num_ancestors_frame.num_ancestors_e[num_ancestors_frame.time .== current_time] .= num_ancestors_e
        num_ancestors_frame.num_ancestors_i[num_ancestors_frame.time .== current_time] .= num_ancestors_i
    end

    return num_ancestors_frame
end

# function find_pair_coal_times(pair_ids, individ_frame)
#     # find the most recent common ancestor of the pair
#     pair1_history = vcat(individ_frame.history[pair_ids[1]], pair_ids[1])
#     pair2_history = vcat(individ_frame.history[pair_ids[2]], pair_ids[2])
#     mrca = maximum(intersect(pair1_history, pair2_history))
#     # two cases, first case, the mrca is one of the pair (i.e. one of the two individuals infected the other)
#     # in this case, the coalescence time is the time of infection of the sampled individual who is not the mrca
#     if mrca == pair_ids[1]
#         mrca_index2 = findfirst(pair2_history .== mrca)
#         coal_time = individ_frame.infec_time[pair2_history[mrca_index2 + 1]]
#     elseif mrca == pair_ids[2]
#         mrca_index1 = findfirst(pair1_history .== mrca)
#         coal_time = individ_frame.infec_time[pair1_history[mrca_index1 + 1]]
#     else 
#         # now the third case is more complicated
#         # in the third case the mrca is not one of the pair members
#         # in which case, the mrca is some kind of ancestor of both pair members, 
#         # for a single pair member, the mrca either infected them, or one of their more direct ancestors
#         # in either case, the coalescence time is the minimum of the infection time of the lineage directly after the mrca
#         # for each pair member 
#         mrca_index1 = findfirst(pair1_history .== mrca)
#         mrca_index2 = findfirst(pair2_history .== mrca)
#         inf_time1 = individ_frame.infec_time[pair1_history[mrca_index1 + 1]]
#         inf_time2 = individ_frame.infec_time[pair2_history[mrca_index2 + 1]]
#         coal_time = minimum([inf_time1, inf_time2])
#     end 
#     return coal_time, mrca
# end

# function find_samp_coal_times(samp_ids, individ_frame)
#     # for all pairs of samp_ids, find the coalescence time
#     active_ids = samp_ids
#     # create vector of length n-1 to store the coalescence times
#     coalescence_times = Vector{Union{Missing, Float64}}(undef, length(samp_ids) -1)
#     for i in 1:(length(samp_ids) -1) 
#     pairs = collect(combinations(active_ids, 2))
#     pair_coal_times = Array{Union{Missing, Float64}}(undef, length(pairs), 2)
#     for (j, pair) in enumerate(pairs)
#         res = find_pair_coal_times(pair, individ_frame)
#         pair_coal_times[j,1] = res[1]
#         pair_coal_times[j,2] = res[2]
#     end
#     # print("active_ids")
#     # print(active_ids)
#     # find the maximum coalescence time (whichever one is sooner happens first)
#     max_coal_time = maximum(pair_coal_times[:,1])
#     # find the pair with the maximum coalescence time
#     max_coal_pair = pairs[argmax(pair_coal_times[:,1])]
#     # find the mrca of the pair with the max coalescence time
#     max_coal_mrca = pair_coal_times[argmax(pair_coal_times[:,1]),2]
#     # remove the pair with the maximum coalescence time from the active_ids
#     active_ids = setdiff(active_ids, max_coal_pair)
#     old_active_ids = copy(active_ids)
#     # now we need to find the ancestors of all of the remaining lineages at the time of the coalescence
#     for old_lineage in old_active_ids
#         # find the vector of possible lineages, which includes the active_id
#         possible_lineages = vcat(individ_frame.history[old_lineage], old_lineage)
#         # filter individ_frame for the possible lineages
#         possible_frame = filter(row -> row.id in possible_lineages, eachrow(individ_frame))
#         # find all lineages who were infected before the coalescence time and who have not yet recovered by the coalescence time
#         eligible_lineages = filter(row -> row.infec_time < max_coal_time && (ismissing(row.recov_time) || row.recov_time > max_coal_time) && row.id in possible_lineages, eachrow(individ_frame))
#         # pick the most recent lineage 
#         closest_lineage = eligible_lineages.id[argmax(eligible_lineages.infec_time)]
#         # remove the active_id from the active_ids
#         active_ids = filter!(x -> x != old_lineage, active_ids)
#         # add the closest lineage to the active_ids
#         push!(active_ids, closest_lineage)
#     end

#     # add the mrca to the active_ids
#     push!(active_ids, max_coal_mrca)
#     # store the minimum coalescence time
#     coalescence_times[i] = max_coal_time
#     end
#     # return the coalescence times
#     return coalescence_times
# end
# """ find the first time when a lineage switches from I to E """
# function find_soonest_i2e_time(I_ids, individ_frame)
#     # find the soonest time when the lineages switch from I to E 
#     #this is simply the time at which they became infectious
#     infectious_times = individ_frame.infectious_time[I_ids]
#     # find the maximum time
#     soonest_time = maximum(infectious_times)
#     # find the id of the individual who became infectious at the soonest time
#     soonest_id = I_ids[argmax(infectious_times)]
#     return soonest_time, soonest_id
# end 
# """ find the first time when a lineage switches from E to I"""
# function find_soonest_e2i_time(E_ids, individ_frame)
#     # find the soonest time when the lineages switch from E to I
#     # this is the time at which they were infected
#     inf_times = individ_frame.infec_time[E_ids]
#     # find the minimum time
#     soonest_time = minimum(inf_times)
#     # find the id of the individual who was infected at the soonest time
#     soonest_id = E_ids[argmin(inf_times)]
#     # find their infector
#     infector = individ_frame.history[soonest_id][end]
#     return soonest_time, soonest_id, infector
# end
# function create_ei_coal_tree(samp_ids, individ_frame, samp_time1)
#     # first find the mrca, the infectious time of the mrca is the time of the root of the tree
#     # intersect all histories of all samp_ids, find the minimum
#     all_histories = [individ_frame.history[samp_id] for samp_id in samp_ids]
#     mrca = minimum(intersect(all_histories...))
#     root_time = individ_frame.infectious_time[mrca]
#     tree = TreeNode(:root, root_time, mrca, "I")
#     childone = TreeNode(:sampled_death, samp_time1, samp_ids[1], "I")
#     tree = attach!(tree, childone)

#     # first find the smallest coalescent time for the samp_ids
#     pairs = collect(combinations(samp_ids, 2))
#     pair_coal_times = Array{Float64}(undef, length(pairs), 2)
#     for (j, pair) in enumerate(pairs)
#         res = find_pair_coal_times(pair, individ_frame)
#         pair_coal_times[j,1] = res[1]
#         pair_coal_times[j,2] = res[2]
#     end
#     # find the soonest coalescence time
#     min_coal_time = maximum(pair_coal_times[:,1])
#     # find the pair with the minimum coalescence time
#     min_coal_pair = pairs[argmax(pair_coal_times[:,1])]
#     # find the mrca of the pair with the minimum coalescence time
#     min_coal_mrca = pair_coal_times[argmax(pair_coal_times[:,1]),2]


    
# end 

# simulate ei model with timevarying alpha
function sim_ei_timevar_alpha(I_init::Int64, alpha_vec, alpha_times, nu::Float64, gamma::Float64, stop_time::Float64)
    # Setting up initial states and frames
    t = 0.0
    num_exposed = 0
    num_infectious = I_init
    num_recovered = 0
    current_alpha = alpha_vec[1]
    rt = current_alpha/nu 
    state_frame = [t num_exposed num_infectious num_recovered current_alpha rt]

    # initialize ids and individ_frame
    exposed_ids = Int[]
    infectious_ids = Int[]
    recovered_ids = Int[]
    # Initialize the last id
    current_id = I_init
    
    # Tracking history and times
    individ_frame = DataFrame(id=Int[], 
                              infec_time=Union{Missing, Float64}[], 
                              infectious_time=Union{Missing, Float64}[], 
                              recov_time=Union{Missing, Float64}[], 
                              history=Any[])  # List of who infected them

    # Initialize initial exposed individuals
    for i in 1:I_init
        push!(infectious_ids, i)
        push!(individ_frame, (i, -1, 0, missing, [i]))
    end
    # simulate till end point reached
    while (num_exposed > 0 || num_infectious > 0) && t < stop_time
        # update alpha
        alpha_idx = if t== 0
            1 
        else 
            findlast(alpha_times .<= t)
        end
        current_alpha = alpha_vec[alpha_idx]
        next_change = if alpha_idx == length(alpha_vec)
            stop_time + 100000
        else 
            alpha_times[alpha_idx + 1]
        end 
        # Time to the next event
        next_event = rand(Exponential(1 / (current_alpha * num_infectious + gamma * num_exposed + nu * num_infectious)))

        # Update time
        # first check if an alpha change has occurred 
        if t + next_event > next_change
            t = next_change
            # update state_frame
            alpha_idx = if t== 0
                1 
            else 
                findlast(alpha_times .<= t)
            end
            current_alpha = alpha_vec[alpha_idx]    
            rt = current_alpha/nu
            state_frame = vcat(state_frame, [t num_exposed num_infectious num_recovered current_alpha rt])
        else
            t += next_event

            # Choose which event happens proportional to the rates
            which_event = sample(["infection", "infectious", "recovery"], Weights([current_alpha * num_infectious, gamma * num_exposed, nu * num_infectious]))

            # If the event is infection
            if which_event == "infection"
                # update current id
                current_id += 1
                # choose an infector from among the infectious individuals at random
                infector_id = sample(infectious_ids)
                # Update states
                num_exposed += 1
                push!(exposed_ids, current_id)
                push!(individ_frame, (current_id, t, missing, missing, unique([individ_frame.history[infector_id]..., infector_id])))
            # If the event is becoming infectious
            elseif which_event == "infectious"
                # choose an exposed individual at random
                exposed_id = sample(exposed_ids)
                # Update states
                num_exposed -= 1
                num_infectious += 1
                push!(infectious_ids, exposed_id)
                # remove exposed_id from exposed_ids
                exposed_ids = filter!(x -> x != exposed_id, exposed_ids)

                individ_frame.infectious_time[exposed_id] = t
            # If the event is recovery
            else
                # choose an infectious individual at random
                infectious_id = sample(infectious_ids)
                # Update states
                num_infectious -= 1
                num_recovered += 1
                push!(recovered_ids, infectious_id)
                # remove infectious_id from infectious_ids
                infectious_ids = filter!(x -> x != infectious_id, infectious_ids)
                individ_frame.recov_time[infectious_id] = t
            end
            # update state_frame
            rt = current_alpha/nu
            state_frame = vcat(state_frame, [t num_exposed num_infectious num_recovered current_alpha rt])
        end 
    end
    return (individ_frame, DataFrame(state_frame, [:time, :E, :I, :R, :alpha, :rt]))
end


