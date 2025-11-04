function sample_init_params(log_rw_mean, log_rw_sigma_sd, log_rt_init_mean, log_rt_init_sd,
    log_gamma_mean, log_gamma_sd, log_nu_mean, log_nu_sd, log_e0_mean, log_e0_sd, log_i0_mean, log_i0_sd, 
  alpha_times, curr_lin, comp_times, reverse_comp_times, est_times)
    E_traj = zeros(length(comp_times))
    I_traj = zeros(length(comp_times))
    rw_sigma = 0.0
    init_rts = zeros(length(alpha_times))
    init_gamma = 0.0
    init_nu = 0.0
    init_e0 = 0.0
    init_i0 = 0.0
    bad_params = true
    total_pop_big_enough = zeros(Bool, length(curr_lin))
    while bad_params
        # sample from the priors 
        rw_sigma = exp(rand(Normal(log_rw_mean, log_rw_sigma_sd)))
        first_rt = exp(rand(Normal(log_rt_init_mean, log_rt_init_sd)))
        diffs = randn(length(alpha_times)-1) .* rw_sigma
        rt_no_init = exp.(log(first_rt) .+ cumsum(vec(diffs)))
        init_rts = vcat(first_rt, rt_no_init)
        init_gamma = exp(rand(Normal(log_gamma_mean, log_gamma_sd)))
        init_nu = exp(rand(Normal(log_nu_mean, log_nu_sd)))
        init_e0 = exp(rand(Normal(log_e0_mean, log_e0_sd)))
        init_i0 = exp(rand(Normal(log_i0_mean, log_i0_sd)))
        alpha_vec = init_rts .* init_nu
        calc_ei_trajectoriesv2!(comp_times, alpha_times, alpha_vec, init_gamma, init_nu, init_e0,init_i0, E_traj, I_traj)
        reverse_E = vcat(reverse(E_traj), init_e0)[2:end]
        reverse_I = vcat(reverse(I_traj), init_i0)[2:end]
        total_pop = reverse_E .+ reverse_I
        # check if the total population size is big enough at all times
        total_pop_big_enough = check_total_pop_size(est_times, curr_lin, reverse_comp_times, reverse_E, reverse_I)
        # but also check if the first value of reverse_I is big enough, 
        # because we know all lineages are initially in I 
        if !any(total_pop .> 8E9) && all(total_pop_big_enough .== 1) && !(reverse_I[1] < curr_lin[1])
            bad_params = false
        end
    end
    return rw_sigma, init_rts, init_gamma, init_nu, init_e0, init_i0
end
"""
check_total_pop_size(total_pop, curr_lin_vec)
this function checks if the total population size is too small at any time point
"""
function check_total_pop_size(est_times, curr_lin, reverse_comp_times, reverse_E, reverse_I)
    est_times_idx = 1
    pop_big_enough = zeros(length(curr_lin))
    for i in 1:length(reverse_comp_times)
        if reverse_comp_times[i] == est_times[est_times_idx]
            if floor(reverse_E[i]) + floor(reverse_I[i]) < curr_lin[est_times_idx]
                pop_big_enough[est_times_idx] = 0.0
            else
                pop_big_enough[est_times_idx] = 1.0
            end
            est_times_idx +=1
        end 
    end
    return pop_big_enough 
end