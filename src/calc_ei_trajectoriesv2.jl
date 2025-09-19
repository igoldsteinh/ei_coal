"""
ei_trajectories(all_times, alpha_times, alpha, gamma, nu, init)
calculate the values of E and I at all_times 
- target_times: vector of times at which to calculate E and I (forward time), assuming the start time is time 0 (not included)
- alpha_times: vector of times at which alpha changes, should include 0 (forward time)
- alpha: vector of alpha values, length pre-specified, alpha/nu = Rt, should include initial alpha at time 0
- gamma: 1/gamma is mean of latent period
- nu: 1/nu is mean of infectious period
- init: initial conditions for E and I, length 2 assuming at time 0
Should return two modified vectors E_traj and I_traj, which are of length target_times
"""
function calc_ei_trajectoriesv2!(target_times, alpha_times, alphas, gamma, nu, e0, i0, E_traj, I_traj)
    params = [gamma, nu, alphas[1], target_times[1]]
    E_traj[1] = g1(params) * e0 + g2(params) * i0
    I_traj[1] = h1(params) * e0 + h2(params) * i0
    new_e0 = E_traj[1]
    new_i0 = I_traj[1]
    alpha_idx = 1
    @inbounds for t in 2:length(target_times)
        # find the right alphas value
        if alpha_idx < length(alpha_times) && alpha_times[alpha_idx + 1] < target_times[t]
            alpha_idx += 1
            params[3] = alphas[alpha_idx]
        end
        # calculate the E and I values at time t
        params[4] = target_times[t] - target_times[t-1]
        E_traj[t] = g1(params) * new_e0 + g2(params) * new_i0
        I_traj[t] = h1(params) * new_e0 + h2(params) * new_i0
        new_e0 = E_traj[t]
        new_i0 = I_traj[t]
    end
end    
