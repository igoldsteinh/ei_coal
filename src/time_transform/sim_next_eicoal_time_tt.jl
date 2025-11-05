"""
sim_next_eicoal_time_tt(t_start; reverse_times, samp_time, comp_traj, gamma, 
    alpha, nE, nI)
simulate the next time in the partially observed EI coalescent using time transformation
assuming that the hazard is piecewise constant 
# Arguments
-t_start: current time 
-reverse_times: vector of times at which the rate changes
-samp_time: last sample time in the sample 
-comp_traj: dataframe of population level E and I values
-gamma: gamma parameter
-alpha: alpha parameter
-nE: number of lineages in state E at time t_start
-nI: number of lineages in state I at time t_start
"""
function sim_next_eicoal_time_tt(t_start; reverse_times, samp_time, comp_traj, gamma, 
    alpha, nE, nI)
    target = rand(Exponential(1))
    # Initialize current time
    t = t_start
    rate = 0.0
    result = 0.0
    # at t_start, we don't know the closes change point, so we find it expensively
    cp_idx = findfirst(reverse_times .> t)
    while target > 0 && sum(reverse_times .> t) > 0
        # Compute the rate for the I2E process
        next_cp = reverse_times[cp_idx]
        eval_point = (t + next_cp)/2
        rate, coal_rate, I2E_rate, E2I_rate, E, I = ei_event_rate_sampler(eval_point; 
            init_eff_frame = comp_traj, alpha = alpha, gamma = gamma, max_time = samp_time, nE = Int(nE), nI = Int(nI))
        if (target - (next_cp - t) * rate) <= 0
            # Calculate time of next event
            result = t + target / (rate)
            target = 0
        else
            # Update target based on how much time has passed with this rate
            target -= (next_cp - t) * rate
            t = next_cp
            # update the index of the next cp, which is predictable, it's simply the next one
            cp_idx += 1
        end 
    end
    # calculate what kind of event it is 
    rate, coal_rate, I2E_rate, E2I_rate, E, I = ei_event_rate_sampler(result; 
    init_eff_frame = comp_traj, alpha = alpha, gamma = gamma, max_time = samp_time, nE = Int(nE), nI = Int(nI))
    event = sample(["I2E", "E2I", "coal"], Weights([I2E_rate, E2I_rate, coal_rate]))
    if result < t_start
        event = "invalid"
    end
    return result, event, E, I  # Return time at which the next event occurs
end