
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
we are now using log_e0 as itself, not log(e0 - 1)
"""
function sample_ess_noalpha(q_cur::Vector{Float64}, l_cur::Float64, cholC::AbstractMatrix, log_prior_means::Vector{Float64}, 
    init_lineages::Int, est_times::Vector{Float64}, coal_times::Vector{Float64}, est_states::Vector{Int}, reverse_samp_times, reverse_samp_lin::AbstractVector, 
    alpha_times::Vector{Float64}, A_matrix::Matrix{Float64}, 
    L_matrix::Matrix{Float64}, my_vec, row_vector, vector_cache, cache_dict, 
    ks_dict, expv_cache_dict, mat_size::Int, E_traj::Vector{Float64}, I_traj::Vector{Float64}, reverse_E::Vector{Float64}, 
    reverse_I::Vector{Float64}, alpha_vec, reverse_alpha_vec::Vector{Float64}, lin_E::Vector{Float64}, 
    lin_I::Vector{Float64}, pop_big_enough::Vector{Float64}, tstep_cutoff)
    if isinf(l_cur)
        print("l_cur")
        print(l_cur)
        error("l_cur is inf")
    end

    # Choose ellipse
    ind_params = cholC * randn(size(cholC)[1])
    # suppose the last of them is the rw parameter
    rw_sigma = exp(ind_params[end-1] + log_prior_means[end-1])
    diffs = randn(length(alpha_times)-1) .* rw_sigma
    # suppose the second to last is the 0 mean log_rt_init parameter
    log_rt_no_init_nomean = ind_params[end] .+ cumsum(vec(diffs))
    v = vcat(ind_params, log_rt_no_init_nomean)
    # Log-likelihood threshold
    u = rand()  # uniform random number in [0, 1]
    logy = l_cur + log(u)
    # println("v")
    # println(v)
    # println("u")
    # println(u)
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
    reverse_alpha_vec .= reverse(alpha_vec)
    calc_ei_trajectoriesv2!(comp_times, alpha_times, alpha_vec, gamma, nu, e0, i0, E_traj, I_traj)
    reverse_E[1:end-1] .= reverse(E_traj)[2:end]
    reverse_E[end] = e0
    reverse_I[1:end-1] .= reverse(I_traj)[2:end]
    reverse_I[end] = i0
    max_pop = maximum(reverse_E .+ reverse_I)

    ll = calc_full_likelihood_noalpha(init_lineages, est_times, est_states,
     reverse_comp_times, coal_times, reverse_samp_times, reverse_samp_lin, gamma, 
    reverse_alpha_vec, reverse_E, reverse_I, max_pop, 
    A_matrix, L_matrix, my_vec, row_vector, vector_cache,
     cache_dict, ks_dict, expv_cache_dict, mat_size, lin_E, lin_I, 
    pop_big_enough, tstep_cutoff)

    iter = 1
    while ll < logy || isnan(ll) 
        # Shrink the bracket and try a new point
        # println("logy")
        # println(logy)
        # println("ll")
        # println(ll)
        # println("check")
        # println(ll < logy)
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
        reverse_alpha_vec .= reverse(alpha_vec)
        calc_ei_trajectoriesv2!(comp_times, alpha_times, alpha_vec, gamma, nu, e0, i0, E_traj, I_traj)
        reverse_E[1:end-1] .= reverse(E_traj)[2:end]
        reverse_E[end] = e0
        reverse_I[1:end-1] .= reverse(I_traj)[2:end]
        reverse_I[end] = i0
        max_pop = maximum(reverse_E .+ reverse_I)

        ll = calc_full_likelihood_noalpha(init_lineages, est_times, est_states,
     reverse_comp_times, coal_times, reverse_samp_times, reverse_samp_lin, gamma, 
    reverse_alpha_vec, reverse_E, reverse_I, max_pop, 
    A_matrix, L_matrix, my_vec, row_vector, vector_cache,
     cache_dict, ks_dict, expv_cache_dict, mat_size, lin_E, lin_I, 
    pop_big_enough, tstep_cutoff)
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
    return q, ll
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