"""
"""
function calc_full_likelihood_noalpha(init_lineages::Int, est_times::AbstractVector{Float64}, est_states::Vector{Int},
     reverse_comp_times::AbstractVector{Float64},
    coal_times::AbstractVector{Float64}, 
     reverse_samp_times::AbstractVector, 
    reverse_samp_lin::AbstractVector, gamma::Float64, 
    reverse_alpha_vec::AbstractVector{Float64}, reverse_E::AbstractVector{Float64}, reverse_I::AbstractVector{Float64}, 
    max_pop, 
    A_matrix::AbstractMatrix{Float64}, L_matrix::AbstractMatrix{Float64}, my_vec::AbstractVector{Float64}, row_vector, vector_cache,
     cache_dict, ks_dict, expv_cache_dict, mat_size::Int, lin_E::AbstractVector{Float64}, lin_I::AbstractVector{Float64}, 
    pop_big_enough::Vector{Float64}, tstep_cutoff)
    # make ll_vec, only for debugging
    # ll_vec = zeros(length(est_times))
        # likelihood checks, pop size too large, or if I is too small
    if max_pop > 8E9
        return -Inf
    end 
    # check if the populations are big enough at each coal and sample time
    pop_big_enough = check_total_pop_size_noalpha!(est_times, lin_E, lin_I, reverse_comp_times, reverse_E, reverse_I, pop_big_enough)
    # if they're ever not, reject immediately
    if any(x -> x == 0, pop_big_enough)
        return -Inf
    end
    # set up active_lineages, active_states, and active_times
    # these are the starting times and states, and the ending times and states for each likelihood chunk
    # they are updated as we move through estimation times
    my_method = ExpMethodHigham2005()
    active_lineages = init_lineages
    active_start_time = 0.0
    active_est_time = est_times[1]
    active_start_state = 1
    active_est_state = est_states[1]
    # set the ticker for sampling vs coalescent times vs alpha times
    coal_time_idx = 1
    next_coal_time = coal_times[coal_time_idx]
    samp_time_idx = 1
    next_samp_time = if length(reverse_samp_times) > 0
        reverse_samp_times[samp_time_idx]
    else 
        -1
    end
    alpha_t = reverse_alpha_vec[1]
    # the first value of reverse_E and reverse_I corresponds to the time of last sampling 
    # so the second value corresponds to the first estimation time
    log_lik = 0.0
    @inbounds for j in 1:(length(est_times) - 1)
        delta_t = active_est_time - active_start_time
        if active_est_time == next_coal_time
            if (active_lineages + 1) > mat_size && delta_t < tstep_cutoff # use the krylov method, will be some error but hopefully not enough to matter
                @views new_ll, alpha_t = calc_coal_ll_noalpha_krylov(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
    L_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], my_vec[1:(active_lineages +1)], ks_dict[active_lineages + 1], 
    expv_cache_dict[active_lineages + 1], active_lineages, active_est_time, active_est_state, active_start_time, active_start_state,
    gamma, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I)
            else # use the matrix method
                @views new_ll, alpha_t = calc_coal_ll_noalpha_higham(active_lineages, active_est_time, active_est_state, 
                active_start_time, active_start_state, gamma, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, 
                reverse_I,  A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], L_matrix[1:(active_lineages + 1), 1:(active_lineages -1)], 
                my_method, cache_dict[(active_lineages + 1)], row_vector[1:(active_lineages + 1)], vector_cache[1:(active_lineages + 1)]) 
            end
            log_lik += new_ll
            # ll_vec[j] = new_ll
            active_lineages -= 1
            coal_time_idx += 1
            next_coal_time = coal_times[coal_time_idx]
        elseif active_est_time == next_samp_time
            if (active_lineages + 1) > mat_size && delta_t < tstep_cutoff # use the krylov method
                @views new_ll, alpha_t = calc_samp_ll_noalpha_krylov(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                ks_dict[active_lineages + 1], expv_cache_dict[active_lineages + 1], active_lineages, active_est_time, active_est_state, 
                active_start_time, active_start_state, gamma, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, 
                reverse_I)
            else # use the matrix method
                @views new_ll, alpha_t = calc_samp_ll_noalpha_higham(active_lineages, active_est_state, active_est_time, 
                active_start_state, active_start_time, gamma, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, 
                reverse_I,  A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], row_vector[1:(active_lineages + 1)], 
                my_method, cache_dict[(active_lineages + 1)], vector_cache[1:(active_lineages + 1)]) 
            end
            log_lik += new_ll
            active_lineages += reverse_samp_lin[samp_time_idx]
            samp_time_idx += 1
            if samp_time_idx <= length(reverse_samp_times)
                next_samp_time = reverse_samp_times[samp_time_idx]
            end
        end
        # update the start and end states
        active_start_state = active_est_state
        active_est_state = est_states[j+1]
        # update the active start time
        active_start_time = active_est_time
        # update the active est time
        active_est_time = est_times[j+1]
    end 
    # now do the the two lineage
    @views  last_prob = calc_coal_twolin_pdf_noalpha(active_est_time, active_start_state, active_start_time, 
    gamma::Float64, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I,  A_matrix[1:3, 1:3], 
     my_method, cache_dict[3], row_vector[1:3],
     vector_cache[1:3])
    log_lik += log(last_prob)
    # ll_vec[end] = log(last_prob)
    return log_lik
end

"""
check_comp_pop_size(lin_E, lin_I, reverse_E, reverse_I)
assuming reverse_E, reverse_I are same length as lin_E, lin_I
check that the populations in the compartments are large enough 
reverse_E and reverse_I are longer than lin_E and lin_I
reverse_comp_times is the same length as reverse_E and reverse_I
est_times is the same length as lin_E and lin_I
"""
function check_total_pop_size_noalpha!(est_times, lin_E, lin_I, reverse_comp_times, reverse_E, reverse_I, pop_big_enough)
    # check that the total population size is large enough
    est_times_idx = 1
    for i in 1:length(reverse_comp_times)
        if reverse_comp_times[i] == est_times[est_times_idx]
            if floor(reverse_E[i]) + floor(reverse_I[i]) < lin_E[est_times_idx] + lin_I[est_times_idx]
                pop_big_enough[est_times_idx] = false
            else
                pop_big_enough[est_times_idx] = true
            end
            est_times_idx +=1
        end 
    end
    return pop_big_enough 
end
"""
"""
function create_currlin(initial_lineages, coal_and_samp_times, coal_times, reverse_samp_lin)
    coal_idx = 1
    samp_idx = 1 
    curr_lin = zeros(length(coal_and_samp_times))
    if coal_and_samp_times[1] == coal_times[coal_idx]
        curr_lin[1] = initial_lineages - 1
        if coal_idx < length(coal_times)
                coal_idx += 1
        end
    else
        curr_lin[1] = initial_lineages + reverse_samp_lin[samp_idx]
        if samp_idx < length(reverse_samp_lin)
                samp_idx += 1
        end
    end 
    for i in 2:length(coal_and_samp_times)
        if coal_and_samp_times[i] == coal_times[coal_idx]
            curr_lin[i] = curr_lin[i-1] - 1
            if coal_idx < length(coal_times)
                coal_idx += 1
            end
    else 
            curr_lin[i] = curr_lin[i-1] + reverse_samp_lin[samp_idx]
            if samp_idx < length(reverse_samp_lin)
                samp_idx += 1
            end
        end
    end 
    return curr_lin 
end 
function create_linElinI(sampled_node_states, curr_lin)
    lin_E = zeros(length(sampled_node_states))
    lin_I = zeros(length(sampled_node_states))
    lin_E .= sampled_node_states .- 1
    lin_I .= curr_lin .- lin_E
    return lin_E, lin_I
end 
