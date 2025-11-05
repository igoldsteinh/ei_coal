"""
propose_ei_coal_tree()
propose a EI coal tree for known samples and known EI counts in forward time 
each new event is simulated via time transformation
if the number of samples exceeds the population size, reject the trajectory and start over
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
function propose_ei_coal_tree_tt_rejectbadpop(max_time; samp_time, reverse_times, comp_traj, gamma, 
    alpha, init_nE, init_nI)
    current_time = 0.0
    nE = init_nE
    nI = init_nI
    coal_times = []
    all_times = [0.0]
    nE_traj = [init_nE]
    nI_traj = [init_nI]
    event_types = ["start"]
    valid = false 
    while  valid == false
        current_time = 0.0
        nE = init_nE
        nI = init_nI
        coal_times = []
        all_times = [0.0]
        nE_traj = [init_nE]
        nI_traj = [init_nI]
        event_types = ["start"]
        while nE + nI > 1 
            # sample the next event time
            current_time, event, E, I = sim_next_eicoal_time_tt(current_time; reverse_times = reverse_times, samp_time = samp_time, 
            comp_traj = comp_traj, gamma = gamma, alpha = alpha, nE = nE, nI = nI)
            if event == "invalid"
                break
            end
            push!(all_times, current_time)
            push!(event_types, event)
            if event == "I2E"
                nE += 1
                nI -= 1
                push!(nE_traj, nE)
                push!(nI_traj, nI)
            elseif event == "E2I"
                nE -= 1
                nI += 1
                push!(nE_traj, nE)
                push!(nI_traj, nI)
            elseif event == "coal"
                nE -= 1
                push!(nE_traj, nE)
                push!(nI_traj, nI)
                push!(coal_times, current_time)
            end
            # check if the populations have gotten too big
            if nE > E || nI > I
                break
            end
            # else if we make it to the end, set valid to true, and get out
            if nE + nI == 1
                valid = true
            end 
        end
    end
    return (
        all_times = all_times,
        event_types = event_types,
        nE_traj = nE_traj,
        nI_traj = nI_traj,
        coal_times = coal_times
    )
end

