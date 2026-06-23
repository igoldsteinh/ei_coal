
using LinearAlgebra
"""
ess(q_cur, l_cur, loglik, cholC; kwargs...)
peform an elliptical slice sampling step
allow for tuning of the initial angle proposal
@params
    q_cur: current state (with 0 mean)
    l_cur: current log-likelihood (full likelihood)
    loglik: function to compute the log-likelihood
    cholC: Cholesky factor of the covariance matrix of the independent parameters
    kwargs: additional arguments for the loglik function
@returns
    q: new state
    u: new log-likelihood
    Ind: indicator of whether the proposal was accepted (1) or rejected (0)
the order of params is log_gamma, log_nu, log_e0, log_i0, log_rw_sigma, log_rt_init,  log rt no init
we are now using log_e0 as itself
USE SKYLINE PRIOR NOT GMRF
"""
function sample_ess_simple_skyline(q_cur::Vector{Float64}, l_cur::Float64, cholC::AbstractMatrix, log_prior_means::Vector{Float64}, 
    num_lineages::Int, est_times::Vector{Float64}, coal_times::Vector{Float64}, est_states::Vector{Int}, start_time::Float64, 
    last_samp_time::Float64, reverse_samp_times::Vector{Float64}, reverse_samp_lin::AbstractVector, 
    alpha_times::Vector{Float64}, ll_vec::Vector{Float64}, lik_vec::Vector{Vector{Float64}}, A_matrix::Matrix{Float64}, 
    L_matrix::Matrix{Float64}, my_vec::Vector{Float64}, my_method, cache_dict, temp::Vector{Float64}, L_vec::Vector{Float64}, 
    ks_dict, expv_cache_dict, mat_size::Int, E_traj::Vector{Float64}, I_traj::Vector{Float64}, reverse_E::Vector{Float64}, 
    reverse_I::Vector{Float64}, total_pop::Vector{Float64}, alpha_vec::Vector{Float64}, lin_E::Vector{Float64}, 
    lin_I::Vector{Float64}, pop_big_enough::Vector{Float64}, tstep_cutoff)
    if isinf(l_cur)
        print("l_cur")
        print(l_cur)
        print("ll_vec")
        print(sum(ll_vec))
        error("l_cur is inf")
    end
    # Choose ellipse
    ind_params = cholC * randn(size(cholC)[1])
    # suppose the last of them is the rw parameter
    rw_sigma = exp(ind_params[end-1] + log_prior_means[end-1])
    diffs = randn(length(alpha_times)-1) .* rw_sigma
    # suppose the second to last is the 0 mean log_rt_init parameter
    # we're doing skyline now, all the rt values have the exact same prior
    # and they're all independent of each other
    log_rt_no_init_nomean = randn(length(alpha_times)-1) .* cholC[6,6]
    v = vcat(ind_params, log_rt_no_init_nomean)
    # Log-likelihood threshold
    u = rand()  # uniform random number in [0, 1]
    logy = l_cur + log(u)
    # Draw an initial proposal, also defining a bracket
    omega = 2 * pi
    init_t = omega * rand()
    t_min = -init_t
    t_max = init_t + omega
    t = rand() * (t_max - t_min) + t_min  # resample t in the current bracket
    # Initial proposal
    q = q_cur * cos(t) + v * sin(t)
    # transform q to the original space 
    # log_gamma, log_nu, log_e0, log_i0, log_rt_init, log_rw_sigma, log rt no init
    gamma = exp(q[1] + log_prior_means[1])
    nu = exp(q[2] + log_prior_means[2])
    e0 = exp(q[3] + log_prior_means[3])
    i0 = exp(q[4] + log_prior_means[4]) 
    rw_sigma = exp(q[5] + log_prior_means[5])
    init_rt = exp(q[6] + log_prior_means[6])
    alpha_vec[1] = init_rt .* nu
    # the other rts, trying to save space 
    alpha_vec[2:end] .= exp.(log_prior_means[6] .+ q[7:end]) .* nu
    calc_ei_trajectoriesv2!(comp_times, alpha_times, alpha_vec, gamma, nu, e0, i0, E_traj, I_traj)
    reverse_E[1:length(E_traj)] .= reverse(E_traj)
    reverse_E[end] = e0
    reverse_I[1:length(I_traj)] .= reverse(I_traj)
    reverse_I[end] = i0
    total_pop .= reverse_E .+ reverse_I

    ll, ll_vec, pdf_vec = calc_node_loglik_simplev2!(num_lineages, est_times, coal_times, est_states, start_time, last_samp_time, 
    reverse_samp_times, reverse_samp_lin, alpha_times, gamma, alpha_vec, reverse_E, reverse_I, total_pop,
    ll_vec, lik_vec, A_matrix, L_matrix, my_vec, my_method, cache_dict, temp, L_vec, ks_dict, expv_cache_dict, mat_size, lin_E, lin_I,
    pop_big_enough, tstep_cutoff)
    iter = 1
    while ll < logy || isnan(ll) 
        # Shrink the bracket and try a new point
        if t < 0
            t_min = t
        else
            t_max = t
        end
        
        t = rand() * (t_max - t_min) + t_min  # resample t in the current bracket
        q = q_cur * cos(t) + v * sin(t)
        # transform q to the original space 
        # log_gamma, log_nu, log_e0, log_i0, log_rt_init, log_rw_sigma, log rt no init
        gamma = exp(q[1] + log_prior_means[1])
        nu = exp(q[2] + log_prior_means[2])
        e0 = exp(q[3] + log_prior_means[3])
        i0 = exp(q[4] + log_prior_means[4]) 
        rw_sigma = exp(q[5] + log_prior_means[5])
        init_rt = exp(q[6] + log_prior_means[6])
        alpha_vec[1] = init_rt .* nu
        # the other rts, trying to save space 
        alpha_vec[2:end] .= exp.(log_prior_means[6] .+ q[7:end]) .* nu
        calc_ei_trajectoriesv2!(comp_times, alpha_times, alpha_vec, gamma, nu, e0, i0, E_traj, I_traj)
        reverse_E[1:length(E_traj)] .= reverse(E_traj)
        reverse_E[end] = e0
        reverse_I[1:length(I_traj)] .= reverse(I_traj)
        reverse_I[end] = i0
        total_pop .= reverse_E .+ reverse_I

        ll, ll_vec, pdf_vec = calc_node_loglik_simplev2!(num_lineages, est_times, coal_times, est_states, start_time, last_samp_time, 
    reverse_samp_times, reverse_samp_lin, alpha_times, gamma, alpha_vec, reverse_E, reverse_I, total_pop,
    ll_vec, lik_vec, A_matrix, L_matrix, my_vec, my_method, cache_dict, temp, L_vec, ks_dict, expv_cache_dict, 
    mat_size, lin_E, lin_I, pop_big_enough, tstep_cutoff)
        iter += 1
    end
    if isinf(ll) || isnan(ll) || ll < logy
        print("ll")
        print(ll)
        print("logy")
        print(logy)
        print("old_ll")
        print(l_cur)
        error("log likelihood is inf or nan or less than logy")
    end
    return q, ll, ll_vec, pdf_vec
end

"""
create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rt_init_sd, log_rw_sigma_sd)
return cholesky decomposition of covariance matrix of independent prior parameters
the order is gamma, nu, e0, i0, rw_sigma, rt_init
"""
# log_rt_init_sd = 0.2
# log_gamma_sd = 0.25
# log_nu_sd = 0.25
# log_e0_sd = 0.05
# log_i0_sd = 0.05
# log_rw_sigma_sd = 0.1
function create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd,  log_rw_sigma_sd,log_rt_init_sd,)
    cov_matrix = Diagonal(vcat(log_gamma_sd^2, log_nu_sd^2, log_e0_sd^2, log_i0_sd^2, log_rw_sigma_sd^2, log_rt_init_sd^2, ))
    cholC = cholesky(cov_matrix)
    return cholC.L
end 

"""
sample_fixed_priors(cholC, log_prior_means, num_samps)
"""
function sample_fixed_priors(cholC, log_prior_means, num_samps)
    # create the matrix of paramater values
    my_samples = zeros(num_samps, length(log_prior_means))
    for i in 1:num_samps
        my_samples[i,:] = cholC * randn(size(cholC)[1]) .+ log_prior_means
    end
    return my_samples
end