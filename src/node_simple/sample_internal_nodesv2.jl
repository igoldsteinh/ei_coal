"""
calc_coal_pdf_simplev2()
calculate the pdf of a single coalescent event for each possible absorbption state
NOTE the resulting pdf is indexed by the end state for a known start state
in contrast in the normal likelihood calculation the pdf is indexed by the start state for a known end state, so the dimensions are different
here we assume that no rate changes occur during the interval
# Arguments
-num_lineages: current number of lineages
-coal_time: time of coalescent
-init_dist: vector, one entry is one which is the initial state
-start_time: time at start of interval
-gamma: gamma param
-alpha_t: value of alpha in the interval
-E: value of E in the interval
-I: value of I in the interval
-A_matrix: contains rates of non-coaelscent states
-L_matrix: contains rates of coalescent
-new_inverse_term: needed matrix for reducing allocations
-pdf: returns pdf vector
-my_method: method of matrix exponentiation
-my_cache: cache used to reduced allocations
-: remove later
"""
function calc_coal_pdf_simplev2(num_lineages::Int, coal_time::Float64, init_dist::AbstractVector{Float64}, start_time::Float64, 
    gamma::Float64, alpha_t::Float64, E::Float64, I::Float64, A_matrix::AbstractMatrix{Float64}, 
     L_matrix::AbstractMatrix{Float64}, new_inverse_term::AbstractMatrix{Float64}, 
    pdf::AbstractVector{Float64}, my_method, my_cache) 
        # initialize the matrices
        new_inverse_term .= 0
        pdf .= 0
        update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, I)
        update_L_matrix!(L_matrix, num_lineages, A_matrix)
        delta_t = coal_time - start_time
        # A_matrix now holds the exponential of the A matrix
        exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
        # this is the product term (exp(At)L
        mul!(new_inverse_term,(A_matrix), L_matrix)
        # this is the product term (exp(At)L)^{T}init_dist
        mul!(pdf, transpose(new_inverse_term), init_dist)
        return pdf
end 
"""
calc_samp_pmf_cachedv2!() 
final vector is indexed by end state for a known start state 
assuming no rate changes during the interval
we may be able to improve this with krylov for larger matrices
calculating the pmf for a sampling or alpha node
# Arguments
-num_lineages: number of lineages
-est_time: time of event
-init_dist: vector, one at initial state
-start_time: time of interval start
-gamma: gamma parameter
-alpha_t: value of alpha in interval
-E: value of E in interval
-I: value of I in interval
-A_matrix: contains rates of non-coalescence
-my_method: method for matrix exponentiation
-my_cache: cache for matrix exponentiation
-temp_vec: needed for saving memory
"""
function calc_samp_pmf_cachedv2!(num_lineages::Int, est_time::Float64, init_dist::AbstractVector{Float64}, start_time::Float64, 
    gamma::Float64, alpha_t::Float64, E::Float64, I::Float64, A_matrix::AbstractMatrix{Float64}, my_method, my_cache, 
    temp_vec::AbstractVector{Float64})
    # find the index of the reverse_time before the est time 
    delta_t = est_time - start_time
    update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, I)
    mul!(temp_vec, transpose(exponential!(rmul!(A_matrix, delta_t), my_method, my_cache)), init_dist)
    return temp_vec
end 
"""
sample_internal_nodesv2!() 
sample the number of lineages in states E and I at est_times
assumes that est_times includes alpha_times, coal_times and samp_times
# Arguments
-sampled_node_states: pre-stored output vector
-num_lineages: number of lineages at time 0 (backwards time)
-est_times: sorted times to estimate states
-coal_times: sorted times of coalescence
-init_dist: vector, 1 at initial state
-start_time: set to 0, re-used
-last_samp_time: time of the last sample in forward time
-reverse_samp_times: sorted reverse sample times, not including last sample time
-reverse_samp_lin: number of lineages sampled at reverse_samp_times
-alpha_times: forward times when alpha changes, includes time 0
-gamma: gamma parameter
-alpha_vec: vector of alphas, corresponding entry is what it turns into at the alpha time
-reverse_E: vector of values of E in reverse time
-reverse_I: vector of values of I in reverse time
-cache_dict: dictionary of caches used in exponentiation
-mat_size: cutoff, above krylov subspaces can be used to save computation time
-ks_dict: dictionary of caches for krylov subspaces
-expv_cache_dict: cache dictionary for expv function
-A_matrix: rate matrix of non-coalescent states
-L_matrix: rate matrix of coalescent states
-new_inverse_term: needed for saving memory in likelihood calcs
-pdf: vector of pdf values
-my_method: method for matrix exponentiation
-temp_vec: used for saving cache
-L_vec: used for calculating pdf with two lineages
-my_vec: used for ll calcs
-ll_vec: vector of likelihod contributions of each est_time
-tstep_cutoff: above this time, krylov subspace will not be used (poor approximation)
"""
function sample_internal_nodesv2!(sampled_node_states::Vector{Int}, num_lineages::Int, est_times::AbstractVector{Float64}, coal_times::AbstractVector{Float64}, 
    init_dist::AbstractVector{Float64}, start_time::Float64, last_samp_time::Float64, reverse_samp_times::AbstractVector{Float64}, 
    reverse_samp_lin::AbstractVector, alpha_times::AbstractVector{Float64}, gamma::Float64, alpha_vec::AbstractVector{Float64}, 
    reverse_E::AbstractVector{Float64}, reverse_I::AbstractVector{Float64}, cache_dict::Dict{Int, Any}, mat_size::Int, ks_dict::Dict{Int, Any}, 
    expv_cache_dict::Dict{Int, Any}, A_matrix::AbstractMatrix{Float64},  L_matrix::AbstractMatrix{Float64}, 
    new_inverse_term::AbstractMatrix{Float64},  pdf::AbstractVector{Float64}, my_method, 
    temp_vec::AbstractVector{Float64}, L_vec::AbstractVector{Float64}, my_vec::AbstractVector{Float64}, 
    ll_vec::AbstractVector{Float64},
    tstep_cutoff)
    # initialize many things
    # initialize the active lineages and start and est times
    active_lineages = num_lineages
    active_start_time = start_time
    active_est_time = est_times[1]
    # reset sampled_node_states
    sampled_node_states .= 0
    sampled_node_states[1] = 1
    sampled_node_states[end] = 1
    # set the ticker for sampling vs coalescent times vs alpha times 
    coal_time_idx = 1
    next_coal_time = coal_times[coal_time_idx]
    samp_time_idx = 1
    next_samp_time = if length(reverse_samp_times) > 0
        reverse_samp_times[samp_time_idx]
    else 
        -1
    end 
    alpha_time_idx = length(alpha_times)
    # this is the next alpha time in reverse time (alpha_times are in forward time and need to be reversed)
    next_alpha_time = abs(alpha_times[alpha_time_idx] - last_samp_time)
    # set the initial values of alpha, E and I
    alpha_t = alpha_vec[alpha_time_idx]
    # the first value of reverse_E and reverse_I corresponds to the time of last sampling 
    # so the second value corresponds to the first estimation time
    comp_idx = 2
    E = reverse_E[comp_idx]
    I = reverse_I[comp_idx]
    log_lik = 0.0
    # set initial distribution to 0 and 1 at the first state
    init_dist .= 0
    init_dist[1] = 1.0
    delta_t = active_est_time - active_start_time
    @views for t in 1:(length(est_times)-1)
        if active_est_time == next_coal_time 
            active_pdf = calc_coal_pdf_simplev2(active_lineages, active_est_time, 
            init_dist[1:(active_lineages + 1)], active_start_time, 
            gamma, alpha_t, E, I, 
            A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
            L_matrix[1:(active_lineages + 1), 1:(active_lineages - 1)], 
            new_inverse_term[1:(active_lineages + 1), 1:(active_lineages - 1)], 
         pdf[1:(active_lineages-1)],
            my_method, cache_dict[active_lineages + 1]) 
            # sample proportional to the pdf 
            state = sample(1:(active_lineages-1), Weights(active_pdf))
            sampled_node_states[t+1] = state
            # we need to check if we're over the krylov size,
            # if we are, we need to recalculate the pdf using krylov, otherwise it won't match 
            if active_lineages + 1 > mat_size && delta_t < tstep_cutoff
                lik_vec = calc_coal_pdf_vecform_krylov_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                L_matrix[1:(active_lineages + 1), 1:(active_lineages - 1)], my_vec[1:(active_lineages + 1)], 
                ks_dict[Int(active_lineages + 1)], expv_cache_dict[Int(active_lineages + 1)],active_lineages, active_est_time, 
                state, active_start_time, gamma, alpha_t, E, I)
                active_start_state = sampled_node_states[t]
                ll_vec[t] = log(lik_vec[active_start_state])
                log_lik += log(lik_vec[active_start_state])
            else 
                # we can just use what we have already calculated
                ll_vec[t] = log(active_pdf[state])
                log_lik += log(active_pdf[state])
            end
            coal_time_idx += 1
            next_coal_time = coal_times[coal_time_idx]
            active_lineages -= 1
        elseif active_est_time == next_samp_time 
            pmf = calc_samp_pmf_cachedv2!(active_lineages, active_est_time, 
            init_dist[1:(active_lineages + 1)], active_start_time, 
            gamma, alpha_t, E, I, A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
            my_method, cache_dict[active_lineages + 1], temp_vec[1:(active_lineages + 1)])
            state = sample(1:(active_lineages + 1), Weights(pmf))
            sampled_node_states[t+1] = state
            if active_lineages + 1 > mat_size && delta_t < tstep_cutoff
                lik_vec = calc_samp_lik_vecform_krylov_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                my_vec[1:(active_lineages + 1)], ks_dict[Int(active_lineages + 1)], expv_cache_dict[Int(active_lineages + 1)], 
                active_lineages, active_est_time, state, active_start_time, gamma, alpha_t, E, I)
                active_start_state = sampled_node_states[t]
                ll_vec[t] = log(lik_vec[active_start_state])
                log_lik += log(lik_vec[active_start_state])
            else 
                ll_vec[t] = log(pmf[state])
                log_lik += log(pmf[state])
            end
            active_lineages += reverse_samp_lin[samp_time_idx]
            samp_time_idx += 1
            if samp_time_idx <= length(reverse_samp_times)
                next_samp_time = reverse_samp_times[samp_time_idx]
            end
        elseif active_est_time == next_alpha_time
            # update the state at the alpha time, then update alpha 
            pmf = calc_samp_pmf_cachedv2!(active_lineages, active_est_time, 
            init_dist[1:(active_lineages + 1)], active_start_time, 
            gamma, alpha_t, E, I, A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
            my_method, cache_dict[active_lineages + 1], temp_vec[1:(active_lineages + 1)])
            state = sample(1:(active_lineages + 1), Weights(pmf))
            sampled_node_states[t+1] = state
            # we need to check if we're over the krylov size,
            # if we are, we need to recalculate the pdf using krylov, otherwise it won't match 
            if active_lineages + 1 > mat_size && delta_t < tstep_cutoff
                lik_vec = calc_samp_lik_vecform_krylov_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                my_vec[1:(active_lineages + 1)], ks_dict[Int(active_lineages + 1)], expv_cache_dict[Int(active_lineages + 1)],
                active_lineages, active_est_time, state, active_start_time, gamma, alpha_t, E, I)
                active_start_state = sampled_node_states[t]
                ll_vec[t] = log(lik_vec[active_start_state])
                log_lik += log(lik_vec[active_start_state])
            else 
                ll_vec[t] = log(pmf[state])
                log_lik += log(pmf[state])
            end
            # update alpha_t
            alpha_time_idx -= 1
            if alpha_time_idx > 0
                next_alpha_time = abs(alpha_times[alpha_time_idx] - last_samp_time)
                alpha_t = alpha_vec[alpha_time_idx]
            else
                next_alpha_time = -1
            end 
        end
        # update the active init dist
        init_dist .=0 
        init_dist[state] = 1.0
        # update the active start time
        active_start_time = active_est_time
        # update the active coal time
        active_est_time = est_times[t + 1]
        # update E and I 
        comp_idx += 1
        E = reverse_E[comp_idx]
        I = reverse_I[comp_idx]
        delta_t = active_est_time - active_start_time
    end
    init_dist = zeros(active_lineages + 1)
    init_dist[Int64(sampled_node_states[end-1])] = 1.0
    last_pdf = calc_twolin_pdf_vecform_simplev2(A_matrix[1:3, 1:3], L_vec, 
    my_vec[1:3], my_method, cache_dict[3], active_est_time, active_start_time, gamma, alpha_t, E, I)
    ll_vec[end] = log(last_pdf[sampled_node_states[end-1]])
    log_lik += log(last_pdf[sampled_node_states[end-1]])
    if isinf(log_lik)
        println("log_lik: ", log_lik)
        println("ll_vec: ", ll_vec)
        println("sampled_node_states: ", sampled_node_states)
        println("reverse_E: ", reverse_E)
        println("reverse_I: ", reverse_I)
        println("alpha_vec: ", alpha_vec)
        println("gamma: ", gamma)
        error("log_lik is inf")
    end
    reverse_alpha_times = reverse(abs.(alpha_times[2:end] .- last_samp_time))
    curr_lin = create_curr_lin_vec(num_lineages, coal_times, reverse_samp_times, reverse_samp_lin, reverse_alpha_times)
    lin_E = zeros(length(curr_lin))
    lin_E = sampled_node_states .- 1
    lin_I = zeros(length(curr_lin))
    lin_I = curr_lin .- lin_E
    return log_lik, ll_vec, sampled_node_states
end 
"""
calc_samp_pmf_krylov_cached! 
final vector is indexed by end state for a known start state 
assuming no rate changes during the interval
we may be able to improve this with krylov for larger matrices
not being used for the moment
"""
# function calc_samp_pmf_krylov_cachedv2!(num_lineages, est_time, init_dist, start_time, 
#     gamma, alpha_t, E, I, A_matrix, my_ks)
#     # find the index of the reverse_time before the est time 
#     delta_t = est_time - start_time
#     update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, I)
#     # calculate action of exponential vector
#     ExponentialUtilities.arnoldi!(my_ks, transpose(A_matrix), init_dist; ishermitian=false)
#     ExponentialUtilities.expv!(init_dist, delta_t, my_ks)
#     return init_dist
# end

"""
DEBUGGING FUNCTION IGNORE

safe version of sample_internal_nodes before we get all crazy
assumes that est_times includes alpha_times, coal_times and samp_times
"""
function sample_internal_nodes_safev2(num_lineages, est_times, coal_times, init_dist, start_time, last_samp_time, 
    reverse_samp_times, reverse_samp_lin, alpha_times, gamma, alpha_vec, reverse_E, reverse_I,
      max_lineages, mat_size)
    # initialize many things 
    A_matrix = zeros( max_lineages+1, max_lineages+1)
    L_matrix = zeros( max_lineages+1, max_lineages-1)
    new_inverse_term = zeros(max_lineages + 1, max_lineages - 1)
    pdf = zeros(max_lineages - 1)
    my_method = ExpMethodHigham2005()
    temp_vec = zeros(max_lineages + 1)
    L_vec = zeros(3)
    my_vec = zeros(max_lineages + 1)
    sizes = 2:(max_lineages + 1)
    ks_dict = preallocate_krylov((mat_size+1):max_lineages + 1)
    expv_cache_dict = preallocate_expv_cache((mat_size+1):max_lineages + 1)
    cache_dict = preallocate_caches(sizes, my_method)

    # initialize the active lineages and start and est times
    active_lineages = num_lineages
    active_start_time = start_time
    active_est_time = est_times[1]
    # reset sampled_node_states
    sampled_node_states = zeros(Int64, length(est_times)+1)
    sampled_node_states[1] = 1
    sampled_node_states[end] = 1
    # set the ticker for sampling vs coalescent times vs alpha times 
    coal_time_idx = 1
    next_coal_time = coal_times[coal_time_idx]
    samp_time_idx = 1
    next_samp_time = if length(reverse_samp_times) > 0
        reverse_samp_times[samp_time_idx]
    else 
        -1
    end 
    alpha_time_idx = length(alpha_times)
    # this is the next alpha time in reverse time (alpha_times are in forward time and need to be reversed)
    next_alpha_time = abs(alpha_times[alpha_time_idx] - last_samp_time)
    # set the initial values of alpha, E and I
    alpha_t = alpha_vec[alpha_time_idx]
    # the first value of reverse_E and reverse_I corresponds to the time of last sampling 
    # so the second value corresponds to the first estimation time
    comp_idx = 2
    E = reverse_E[comp_idx]
    I = reverse_I[comp_idx]

    log_lik = 0.0
    # set initial distribution to 0 and 1 at the first state
    init_dist .= 0
    init_dist[1] = 1.0
    # just for troubleshooting
    log_lik_vec = zeros(length(est_times))
    delta_t = active_est_time - active_start_time
    @views for t in 1:(length(est_times)-1)
        if active_est_time == next_coal_time 
            active_pdf = calc_coal_pdf_simplev2(active_lineages, active_est_time, 
            init_dist[1:(active_lineages + 1)], active_start_time, 
            gamma, alpha_t, E, I, 
            A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
            L_matrix[1:(active_lineages + 1), 1:(active_lineages - 1)], 
            new_inverse_term[1:(active_lineages + 1), 1:(active_lineages - 1)], 
             pdf[1:(active_lineages-1)],
            my_method, cache_dict[active_lineages + 1]) 
            # sample proportional to the pdf 
            state = sample(1:(active_lineages-1), Weights(active_pdf))
            sampled_node_states[t+1] = state
            # we need to check if we're over the krylov size,
            # if we are, we need to recalculate the pdf using krylov, otherwise it won't match 
            if active_lineages + 1 > mat_size && delta_t < 1.0  
                lik_vec = calc_coal_pdf_vecform_krylov_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                L_matrix[1:(active_lineages + 1), 1:(active_lineages - 1)], my_vec[1:(active_lineages + 1)], 
                ks_dict[Int(active_lineages + 1)], expv_cache_dict[Int(active_lineages + 1)], active_lineages, active_est_time, 
                state, active_start_time, gamma, alpha_t, E, I)
                active_start_state = sampled_node_states[t]
                log_lik_vec[t] = log(lik_vec[active_start_state])
                log_lik += log(lik_vec[active_start_state])
            else 
                # we can just use what we have already calculated
                log_lik_vec[t] = log(active_pdf[state])
                log_lik += log(active_pdf[state])
            end
            coal_time_idx += 1
            next_coal_time = coal_times[coal_time_idx]
            active_lineages -= 1
        elseif active_est_time == next_samp_time 
            pmf = calc_samp_pmf_cachedv2!(active_lineages, active_est_time, 
            init_dist[1:(active_lineages + 1)], active_start_time, 
            gamma, alpha_t, E, I, A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
            my_method, cache_dict[active_lineages + 1], temp_vec[1:(active_lineages + 1)])
            state = sample(1:(active_lineages + 1), Weights(pmf))
            sampled_node_states[t+1] = state
            if active_lineages + 1 > mat_size && delta_t < 1.0
                lik_vec = calc_samp_lik_vecform_krylov_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                my_vec[1:(active_lineages + 1)], ks_dict[Int(active_lineages + 1)], expv_cache_dict[Int(active_lineages + 1)], active_lineages, 
                active_est_time, state, active_start_time, gamma, alpha_t, E, I)
                active_start_state = sampled_node_states[t]
                log_lik_vec[t] = log(lik_vec[active_start_state])
                log_lik += log(lik_vec[active_start_state])
            else 
                log_lik_vec[t] = log(pmf[state])
                log_lik += log(pmf[state])
            end
            active_lineages += reverse_samp_lin[samp_time_idx]
            samp_time_idx += 1
            if samp_time_idx <= length(reverse_samp_times)
                next_samp_time = reverse_samp_times[samp_time_idx]
            end
        elseif active_est_time == next_alpha_time
            # update the state at the alpha time, then update alpha 
            pmf = calc_samp_pmf_cachedv2!(active_lineages, active_est_time, 
            init_dist[1:(active_lineages + 1)], active_start_time, 
            gamma, alpha_t, E, I, A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
            my_method, cache_dict[active_lineages + 1], temp_vec[1:(active_lineages + 1)])
            state = sample(1:(active_lineages + 1), Weights(pmf))
            sampled_node_states[t+1] = state
            # we need to check if we're over the krylov size,
            # if we are, we need to recalculate the pdf using krylov, otherwise it won't match 
            if active_lineages + 1 > mat_size && delta_t < 1.0
                lik_vec = calc_samp_lik_vecform_krylov_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                my_vec[1:(active_lineages + 1)], ks_dict[Int(active_lineages + 1)], expv_cache_dict[Int(active_lineages + 1)], active_lineages, 
                active_est_time, state, active_start_time, gamma, alpha_t, E, I)
                active_start_state = sampled_node_states[t]
                log_lik_vec[t] = log(lik_vec[active_start_state])
                log_lik += log(lik_vec[active_start_state])
            else 
                log_lik_vec[t] = log(pmf[state])
                log_lik += log(pmf[state])
            end
            # update alpha_t
            alpha_time_idx -= 1
            if alpha_time_idx > 0
                next_alpha_time = abs(alpha_times[alpha_time_idx] - last_samp_time)
                alpha_t = alpha_vec[alpha_time_idx]
            else
                next_alpha_time = -1
            end 
        end
        # update the active init dist
        init_dist .=0 
        init_dist[state] = 1.0
        # update the active start time
        active_start_time = active_est_time
        # update the active coal time
        active_est_time = est_times[t + 1]
        # update E and I 
        comp_idx += 1
        E = reverse_E[comp_idx]
        I = reverse_I[comp_idx]
        delta_t = active_est_time - active_start_time
    end
    init_dist = zeros(active_lineages + 1)
    init_dist[Int64(sampled_node_states[end-1])] = 1.0
    last_pdf = calc_twolin_pdf_vecform_simplev2(A_matrix[1:3, 1:3], L_vec, 
    my_vec[1:3], my_method, cache_dict[3], active_est_time, active_start_time, gamma, alpha_t, E, I)
    log_lik_vec[end] = log(last_pdf[sampled_node_states[end-1]])
    log_lik += log(last_pdf[sampled_node_states[end-1]])
    if isinf(log_lik)
        println("log_lik: ", log_lik)
        println("log_lik_vec: ", log_lik_vec)
        println("sampled_node_states: ", sampled_node_states)
        println("reverse_E: ", reverse_E)
        println("reverse_I: ", reverse_I)
        error("log_lik is inf")
    end
    return log_lik, log_lik_vec, sampled_node_states
end 

