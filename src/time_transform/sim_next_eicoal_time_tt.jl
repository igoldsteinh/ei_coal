"""
function for computing next event times for death processes using using time transformation
assuming that the hazard is piecewise constant 
t_start is the start time 
target is a random variable with distribution exponentiontal(1)
all other variables are needed for the rate functions
model_type is either "I2E" or "E2I" to choose which point process rate function to use
"""
function sim_next_eicoal_time_tt(t_start; reverse_times, samp_time, comp_traj, gamma, 
    alpha, nE, nI)
    target = rand(Exponential(1))
    original_target = target  # Store the original target for debugging
    # Initialize current time
    t = t_start
    rate = 0.0
    result = 0.0
    # create vector of change points which is the merged vectors of 
    # comp_traj.reverse_time, E2I_reverse_time, I2E_reverse_time, coal_times
    # and alpha_times
    # at t_start, we don't know the closes change point, so we find it expensively
    cp_idx = findfirst(reverse_times .> t)
    while target > 0 && sum(reverse_times .> t) > 0
        # Compute the rate for the I2E process
        next_cp = reverse_times[cp_idx]
        eval_point = (t + next_cp)/2
        rate, coal_rate, I2E_rate, E2I_rate = ei_event_rate_sampler(eval_point; 
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
    # if result < t_start
    #     println("t_start: ", t_start, " result: ", result)
    #     println("nE: ", nE, " nI: ", nI)
    #     println("original_target: ", original_target)
    #     error("Result time is less than the start time. This should not happen.")
    # end
    # calculate what kind of event it is 
    rate, coal_rate, I2E_rate, E2I_rate = ei_event_rate_sampler(result; 
    init_eff_frame = comp_traj, alpha = alpha, gamma = gamma, max_time = samp_time, nE = Int(nE), nI = Int(nI))
    event = sample(["I2E", "E2I", "coal"], Weights([I2E_rate, E2I_rate, coal_rate]))
    if result < t_start
        event = "invalid"
    end
    return result, event  # Return time at which the next event occurs
end