"""
sample_nodescdf_andparams!()
for a fixed phylogeny, jointly sample the number of lineages in E and I at internal nodes, sample times, 
    and alpha change points, and also the EI population parameters
# Arguments
-q_cur: log scale mean 0 parameters for ESS sampling
-l_cur: current value of the log-likelihood for ESS sampling
-num_lineages: number of lineages at backwards time 0
-est_times: times when lineage states must be estimated (reverse time)
-coal_times: times of coalescence (reverse time)
-est_states: initial lineage states at est_times
-start_time: set to 0, used internally
-last_samp_time: time of the last sample (forward time)
-reverse_samp_times: sorted sample times in reverse time, excluding time 0
-reverse_samp_lin: number of lineages sampled at reverse_samp_times
-alpha_times: vector of times in forward time in which alpha changes, includes time 0
-mat_size: cutoff at which point Krylov subspace methods are used when possible, default 50
-curr_lin: vector of total lineages at est times
-num_samples: number of total MCMC samples 
-discard_init: number of samples to discard as burn-in
-num_thin: every num_thinth sample is kept 
-tstep_cutoff: above this time, Krylov subspace is not used regardless of mat_size,
1.0 seems to work ok for small samples, 0.5 safer for >=100
"""
function sample_nodescdf_andparams!(q_cur::Vector{Float64}, l_cur::Float64, cholC::AbstractMatrix, log_prior_means::Vector{Float64}, 
    num_lineages::Int, est_times::Vector{Float64}, coal_times::Vector{Float64}, est_states::Vector{Int}, start_time::Float64, 
    last_samp_time::Float64, reverse_samp_times::Vector{Float64}, reverse_samp_lin::AbstractVector, 
    alpha_times::Vector{Float64}, mat_size::Int, curr_lin::Vector{Float64}, num_samples::Int, discard_initial::Int, num_thin::Int,
    tstep_cutoff)
     # idiot proofing 
     if any(curr_lin .<= 0)
         error("curr_lin has negative or 0 values, this is not allowed")
     end
    # create the matrix of paramater values
    my_samples = zeros(Int((num_samples - discard_initial)/num_thin), length(q_cur) + 2)
    # create the matrix of states
    my_states = zeros(Int((num_samples - discard_initial)/num_thin), length(est_states) + 2)
    # create needed matrices and vectors to reduce memory usage 
    ll_vec = Vector{Float64}(undef, length(est_times))
    lik_vec = Vector{Vector{Float64}}(undef, length(est_times))
    max_lineages = Int(maximum(curr_lin))
    sizes = 2:(max_lineages + 1)
    my_method = ExpMethodHigham2005()
    cache_dict = preallocate_caches(sizes, my_method)
    ks_dict = preallocate_krylov(mat_size:max_lineages + 1)
    expv_cache_dict = preallocate_expv_cache((mat_size+1):max_lineages + 1)
    A_matrix = zeros( max_lineages+1, max_lineages+1)
    L_matrix = zeros( max_lineages+1, max_lineages-1)
    new_inverse_term = zeros(max_lineages + 1, max_lineages - 1)
    pdf = zeros(max_lineages - 1)
    my_method = ExpMethodHigham2005()
    temp_vec = zeros(max_lineages + 1)
    L_vec = zeros(3)
    my_vec = zeros(max_lineages + 1)
    temp = similar(my_vec)
    pop_big_enough = zeros(length(curr_lin))


    # create other variable vectors
    alpha_vec = zeros(length(alpha_times))
    reverse_E = zeros(length(comp_times) + 1)
    reverse_I = zeros(length(comp_times) + 1)
    E_traj = zeros(length(comp_times))
    I_traj = zeros(length(comp_times))
    total_pop = zeros(length(comp_times) + 1)
    rt_no_init = zeros(length(alpha_times) - 1) 
    gamma = exp(q_cur[1] + log_prior_means[1])
    nu = exp(q_cur[2] + log_prior_means[2])
    e0 = exp(q_cur[3] + log_prior_means[3])
    i0 = exp(q_cur[4] + log_prior_means[4]) 
    init_rt = exp(q_cur[6] + log_prior_means[6])
    rt_no_init .= exp.(log_prior_means[6] .+ q_cur[7:end])
    alpha_init = init_rt .* nu
    alpha_vec[1] = alpha_init
    alpha_vec[2:end] .= rt_no_init .* nu
    calc_ei_trajectoriesv2!(comp_times, alpha_times, alpha_vec, gamma, nu, e0, i0, E_traj, I_traj)
    reverse_E[1:length(E_traj)] .= reverse(E_traj)
    reverse_E[end] = e0
    reverse_I[1:length(I_traj)] .= reverse(I_traj)
    reverse_I[end] = i0
    total_pop .= reverse_E .+ reverse_I
    lin_E = zeros(length(est_states))
    lin_I = zeros(length(est_states))
    j = 1
    m = 1
    for i in 1:num_samples 
        # sample the parameters
        if i % 1000 == 0
            println("i: ", i, " of ", num_samples)
        end
        # update lineages 
        lin_E .= est_states .- 1
        lin_I .= curr_lin .- lin_E
        q_cur, l_cur, ll_vec, lik_vec = sample_ess_simple(q_cur, l_cur, cholC, log_prior_means, num_lineages, est_times, coal_times, 
        est_states, start_time, last_samp_time, reverse_samp_times, reverse_samp_lin,
        alpha_times, ll_vec, lik_vec, A_matrix, L_matrix, my_vec, my_method, cache_dict, temp, L_vec, ks_dict, expv_cache_dict, 
        mat_size, E_traj, I_traj, reverse_E, reverse_I, total_pop, alpha_vec, lin_E, lin_I, pop_big_enough, tstep_cutoff)
        gamma = exp(q_cur[1] + log_prior_means[1])
        nu = exp(q_cur[2] + log_prior_means[2])
        e0 = exp(q_cur[3] + log_prior_means[3])
        i0 = exp(q_cur[4] + log_prior_means[4]) 
        rw_sigma = exp(q_cur[5] + log_prior_means[5])
        init_rt = exp(q_cur[6] + log_prior_means[6])
        rt_no_init .= exp.(log_prior_means[6] .+ q_cur[7:end])
        if i > discard_initial && i % num_thin == 0
            my_samples[j,1:end-2] .= vcat(gamma, nu, e0, i0, rw_sigma, init_rt, rt_no_init)
            my_samples[j,end-1] = l_cur
            my_samples[j,end] = i
            j += 1
        end
        # update gamma and alphae_vec
        alpha_init = init_rt * nu
        alpha_vec[1] = alpha_init
        alpha_vec[2:end] .= rt_no_init .* nu
        calc_ei_trajectoriesv2!(comp_times, alpha_times, alpha_vec, gamma, nu, e0, i0, E_traj, I_traj)
        reverse_E[1:length(E_traj)] .= reverse(E_traj)
        reverse_E[end] = e0
        reverse_I[1:length(I_traj)] .= reverse(I_traj)
        reverse_I[end] = i0
        total_pop .= reverse_E .+ reverse_I
        # sample all states
        
        l_cur, ll_vec, est_states = sample_internal_nodesv2!(num_lineages, est_times, coal_times, init_dist, 
        start_time, last_samp_time, reverse_samp_times, reverse_samp_lin, alpha_times, gamma, alpha_vec, reverse_E, reverse_I,
        cache_dict, mat_size, ks_dict, expv_cache_dict, A_matrix, L_matrix, 
        new_inverse_term, pdf, my_method, temp_vec, L_vec, my_vec, est_states, ll_vec, tstep_cutoff)
        # store the samples
        if i > discard_initial && i % num_thin == 0
            my_states[m,1:end-2] .= est_states
            my_states[m,end-1] = l_cur
            my_states[m,end] = i
            m += 1
        end
    end
    return my_samples, my_states
end


