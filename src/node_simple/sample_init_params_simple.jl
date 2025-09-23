# function for finding a legal set of parameters
# given a set of initial states 
"""
sample_initial_params_simple()
sample initial parameter values from the prior in such a way that the likelihood won't be zero
(total pop size smaller than earth's, total lineages larger than total pop size)
returns natural scale initial parameter values
# Arguments
-log_rw_mean: mean parameter of log-normal rw prior
-log_rw_sigma_sd: sd parameter of log-normal rw prior
-log_rt_init_mean: mean paramater of log-normal initial rt prior
-log_rt_init_sd: sd parmater of log-normal initial rt prior
-log_gamma_mean: mean paramater of log-normal gamma prior
-log_gamma_sd: sd parameter of log-normal gamma prior
-log_nu_mean: mean parameter of log-normal nu prior
-log_nu_sd: sd parameter of log-normal nu prior
-log_e0_mean: mean parameter of log-normal initial E prior
-log_e0_sd: sd parameter of log-normal intiail E prior
-log_i0_mean: mean parameter of log-normal intial I prior
-log_i0_sd: sd parameter of log-normal intial I prior
-alpha_times: vector forward times of alpha change points
-comp_times: times for which E and I values are known
-curr_lin: total number of lineages at comp times
"""
function sample_initial_params_simple(log_rw_mean, log_rw_sigma_sd, log_rt_init_mean, log_rt_init_sd,
    log_gamma_mean, log_gamma_sd, log_nu_mean, log_nu_sd, log_e0_mean, log_e0_sd, log_i0_mean, log_i0_sd, 
  alpha_times, comp_times, curr_lin)
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
        reverse_E = vcat(reverse(E_traj), init_e0)
        reverse_I = vcat(reverse(I_traj), init_i0)
        total_pop = reverse_E .+ reverse_I
        # check if the total population size is big enough at all times
        total_pop_big_enough = check_total_pop_size(reverse_E, reverse_I, curr_lin)
        # but also check if the first value of reverse_I is big enough, 
        # because we know all lineages are initially in I 
        if !any(total_pop .> 8E9) && all(total_pop_big_enough .== 1) && !(reverse_I[1] < curr_lin[1])
            bad_params = false
        end
    end
    return rw_sigma, init_rts, init_gamma, init_nu, init_e0, init_i0
end
"""
create_curr_lin_vec
this is a function for creating a vector of the number of lineages at each time point
for which the likelihood is calculated (coal, samp, and alpha times)
AFTER THE EVENT HAS OCCURRED
so if the event is a coal event, the number of lineages should go down by 1
this vector is used as part of a likelihood calculation check to ensure legal trajectories quickly
reverse_alpha_times should not include the first alpha time which is always 0
the first value of curr_lin_vec is the number of lineages at the start time
# Arguments
-starting_lineages: number of lineages at backwards time 0
-coal_times: coalescence times in backwards time
-reverse_samp_times: reverse sample times (not including time 0)
-reverse_alpha_times: reverse times that alpha changes
"""
function create_curr_lin_vec(starting_lineages, coal_times, reverse_samp_times, reverse_samp_lin, reverse_alpha_times)
    curr_lin_vec = zeros(length(coal_times) + length(reverse_samp_times) + length(reverse_alpha_times) + 1)
    curr_lin_vec[1] = starting_lineages
    coal_idx = 1
    samp_idx = 1 
    alpha_idx = 1
    next_coal_time = coal_times[coal_idx]
    next_samp_time = if length(reverse_samp_times) > 0
        reverse_samp_times[samp_idx]
    else 
        -1
    end
    next_alpha_time = reverse_alpha_times[alpha_idx]
    # iterate through all times, and update the number of lineages based on the event type
    all_times = sort(vcat(coal_times, reverse_samp_times, reverse_alpha_times))
    for i in 1:length(all_times)
        if all_times[i] == next_coal_time
            curr_lin_vec[i + 1] = curr_lin_vec[i] - 1
            if coal_idx < length(coal_times)
                coal_idx += 1
                next_coal_time = coal_times[coal_idx]
            end
        elseif all_times[i] == next_samp_time
            curr_lin_vec[i + 1] = curr_lin_vec[i] + reverse_samp_lin[samp_idx]
            if samp_idx < length(reverse_samp_lin)
                samp_idx += 1
                next_samp_time = reverse_samp_times[samp_idx]
            end
        elseif all_times[i] == next_alpha_time
            curr_lin_vec[i + 1] = curr_lin_vec[i] 
            if alpha_idx < length(reverse_alpha_times)
                alpha_idx += 1
                next_alpha_time = reverse_alpha_times[alpha_idx]
            end
        end
    end 
    return curr_lin_vec
end
"""
check_total_pop_size(total_pop, curr_lin_vec)
for now we will assume that total_pop and curr_lin_vec are the same length
this function checks if the total population size is too small at any time point
# Arguments
-reverse_E: vector of E values in reverse order
-reverse_I: vector of I values in reverse order
-curr_lin_vec: vector of total lineages 
"""
function check_total_pop_size(reverse_E, reverse_I, curr_lin_vec)
    check_vec = zeros(length(reverse_E))
    for i in 1:length(reverse_E)
        if floor(reverse_E[i]) + floor(reverse_I[i]) < curr_lin_vec[i]
            check_vec[i] = 0
        else
            check_vec[i] = 1
        end
    end
    return check_vec
end