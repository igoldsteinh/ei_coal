
using LinearAlgebra
"""
sample_ess_simple()
peform an elliptical slice sampling step
returns new uncentered parameters, log likelihood, vector of ll contributions, vector of vector of pdfs
# Arguments
-q_cur: current state (with 0 mean)
-l_cur: current log-likelihood (full likelihood)
-cholC: Cholesky factor of the covariance matrix of the independent parameters
-log_prior_means: vector of prior means on log scale
-num_lineages: number of starting lineages
-est_times: times that contribute to the log likelihood
-coal_times: coalescence times
-est_states: number of lineages in states I and E at the est_times
-start_time: zero
-last_samp_time: time of the last sample
-reverse_samp_times: sorted reverse samp times not including the last one
-reverse_samp_lin: number of lineags sampled at reverse_samp_times
-alpha_times: forward times alpha changes
-ll_vec: vector of likelihood contributions of est_times
-lik_vec: vector of likelihood contribution used for cache savings
-A_matrix: rates of non-coalescent states
-L_matrix: rates of coalescent states
-my_vec: used for storing values
-my_method: method of matrix exponenitation
-cache_dict: dictionary of caches for matrix exponentiation
-temp: used for storing values
-L_vec: used for when there are two lineages only
-ks_dict: dictionary of caches for krylov subspace
-expv_cache_dict: dictionary of caches for expv
-mat_size: cutoff above which krylov subpsace can be used
-E_traj: vector of E values
-I_traj: vector of I values
-reverse_E: vector of E values reversed
-reverse_I: vector of I values reversed
-total_pop: sum of E and I
-alpha_vec: vector of alpha values
-lin_E: vector of number of E lineages at est_times
-lin_I: same as lin_E but for I lineages
-pop_big_enough: vector of booleans checking if there are more population members than lineage members
-tstep_cutoff: if time diff is larger than this, krylov subspace is not allowed
the order of params is log_gamma, log_nu, log_e0, log_i0, log_rw_sigma, log_rt_init,  log rt no init
we are now using log_e0 as itself, not log(e0 - 1)
"""
function sample_ess_simple(q_cur::Vector{Float64}, l_cur::Float64, cholC::AbstractMatrix, log_prior_means::Vector{Float64}, 
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
    log_rt_no_init_nomean = ind_params[end] .+ cumsum(vec(diffs))
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
# Arguments
-log_gamma_mean: mean paramater of log-normal gamma prior
-log_gamma_sd: sd parameter of log-normal gamma prior
-log_nu_mean: mean parameter of log-normal nu prior
-log_nu_sd: sd parameter of log-normal nu prior
-log_e0_mean: mean parameter of log-normal initial E prior
-log_e0_sd: sd parameter of log-normal intiail E prior
-log_i0_mean: mean parameter of log-normal intial I prior
-log_i0_sd: sd parameter of log-normal intial I prior
-log_rw_mean: mean parameter of log-normal rw prior
-log_rw_sigma_sd: sd parameter of log-normal rw prior
-log_rt_init_mean: mean paramater of log-normal initial rt prior
-log_rt_init_sd: sd parmater of log-normal initial rt prior
"""
function create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd,  log_rw_sigma_sd,log_rt_init_sd,)
    cov_matrix = Diagonal(vcat(log_gamma_sd^2, log_nu_sd^2, log_e0_sd^2, log_i0_sd^2, log_rw_sigma_sd^2, log_rt_init_sd^2, ))
    cholC = cholesky(cov_matrix)
    return cholC.L
end 

"""
sample_fixed_priors(cholC, log_prior_means, num_samps)
sample the fixed priors
# Arguments
-cholC: cholesky matrix of the variance matrix of the fixed params
-log_prior_means: vector of prior means on log scale
-num_samps: how many samples of the prior?
"""
function sample_fixed_priors(cholC, log_prior_means, num_samps)
    # create the matrix of paramater values
    my_samples = zeros(num_samps, length(log_prior_means))
    for i in 1:num_samps
        my_samples[i,:] = cholC * randn(size(cholC)[1]) .+ log_prior_means
    end
    return my_samples
end