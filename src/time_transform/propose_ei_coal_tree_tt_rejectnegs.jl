"""
propose_ei_coal_tree()
propose a EI coal tree for known samples and known EI counts in forward time 
each new event is simulated via time transformation
if event times are negative, reject and return to the last accepted coalescent time
# Arguments:
-max_time: the time to stop at, not actually necessary
-samp_time: time of the last sample
-reverse_times: vector of times rates change
-comp_traj: frame with E and I counts
-gamma: gamma parameter
-alpha: alpha parameter
-init_nE: initial lineages in state E
-init_nI: initial lineages in state I
"""
function propose_ei_coal_tree_tt_rejectnegs(max_time; samp_time, reverse_times, comp_traj, gamma, 
    alpha, init_nE, init_nI)
    current_time = 0.0
    nE = init_nE
    nI = init_nI
    coal_times = []
    all_times = [0.0]
    nE_traj = [init_nE]
    nI_traj = [init_nI]
    event_types = ["start"]
    while  nE + nI > 1
        potential_events = []
        potential_times = []
        potential_nE = []
        potential_nI = []
        valid = false 
        current_time = all_times[end]
        nE = nE_traj[end]
        nI = nI_traj[end]
        while valid == false 
            # sample the next event time
            current_time, event = sim_next_eicoal_time_tt(current_time; reverse_times = reverse_times, samp_time = samp_time, 
            comp_traj = comp_traj, gamma = gamma, alpha = alpha, nE = nE, nI = nI)
            if event == "invalid"
                break
            end
            push!(potential_times, current_time)
            push!(potential_events, event)
            if event == "I2E"
                nE += 1
                nI -= 1
                push!(potential_nE, nE)
                push!(potential_nI, nI)
            elseif event == "E2I"
                nE -= 1
                nI += 1
                push!(potential_nE, nE)
                push!(potential_nI, nI)
            elseif event == "coal"
                nE -= 1
                push!(potential_nE, nE)
                push!(potential_nI, nI)
                # we reached a coal time, bank all of the potentials and set valid to true
                valid = true
                push!(coal_times, current_time)
                all_times = vcat(all_times, potential_times)
                event_types = vcat(event_types, potential_events)
                nE_traj = vcat(nE_traj, potential_nE)
                nI_traj = vcat(nI_traj, potential_nI)
            end
        end
    end
    # println("coal_times: ", coal_times, "\n")
    return (
        all_times = all_times,
        event_types = event_types,
        nE_traj = nE_traj,
        nI_traj = nI_traj,
        coal_times = coal_times
    )
end

