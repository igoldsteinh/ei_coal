"""
calc_coal_pdf_noalpha
calculate the pdf of a single coalescent event for each possible absorbption state
any changes in alpha are absorbed into the probability
the only possible end states are coalescent times and sampling times
"""
function calc_coal_pdf_noalpha(num_lineages::Int, end_time::Float64, start_state, start_time::Float64, 
    gamma::Float64, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I,  A_matrix::AbstractMatrix{Float64}, 
     L_matrix::AbstractMatrix{Float64}, my_output::AbstractVector{Float64}, my_method, my_cache, row_vector::AbstractVector{Float64},
     vector_cache) 
    # initialize the vectors
    my_output .= 0
    row_vector .= 0
    vector_cache .= 0
    # find the index of the reverse_time before the est time 
    is_between = (reverse_comp_times .> start_time) .& (reverse_comp_times .<= end_time)
    subset_times = reverse_comp_times[is_between]
    subset_E = reverse_E[is_between]
    subset_I = reverse_I[is_between]
    # check if subset_times is of length(1), in which case go the simple route
    # there are no rate changes in between the start and end of the interval
    if length(subset_times) == 1
        E = subset_E[1]
        compI = subset_I[1]
        update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, compI)
        update_L_matrix!(L_matrix, num_lineages, A_matrix)
        delta_t = end_time - start_time
        # A_matrix now holds the exponential of the A matrix
        exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
        # this is the product term (exp(At)L
        mul!(my_output, transpose(L_matrix), A_matrix[start_state,:])
    elseif length(subset_times) > 1
        # now we need to do something more complicated
        # there is at least one alpha update in here
        # go in forward time b/c it means less matrix matrix
        comp_idx = 1
        E = subset_E[comp_idx]
        compI = subset_I[comp_idx]
        update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, compI)
        # update_L_matrix!(L_matrix, num_lineages, A_matrix)
        delta_t = subset_times[1] - start_time
        exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
        row_vector .= @view A_matrix[start_state,:]
        # this first event must have been alpha time, so update alpha_idx
        alpha_idx = findfirst(isequal(alpha_t), reverse_alpha_vec)
        alpha_idx += 1
        for t in 2:length(subset_times)
            alpha_t = reverse_alpha_vec[alpha_idx]
            E = subset_E[t]
            compI = subset_I[t]
            update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, compI)
            delta_t = subset_times[t] - subset_times[t-1]
            if t == length(subset_times)
                # this is a coal_time
                update_L_matrix!(L_matrix, num_lineages, A_matrix)
                exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
                #tranpose(expAt) * row_vector
                mul!(vector_cache, transpose(A_matrix), row_vector)
                # transpose(L) * row_vector
                mul!(my_output, transpose(L_matrix), vector_cache)
            elseif t < length(subset_times)
                exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
                mul!(vector_cache, transpose(A_matrix), row_vector)
                row_vector .= vector_cache
                # this was another alpha event tick up alpha
                alpha_idx += 1
            end 
        end 
    end 
    return my_output, alpha_t
end 
"""
calc_coal_twolin_pdf_noalpha
"""
function calc_coal_twolin_pdf_noalpha(coal_time::Float64, start_state, start_time::Float64, 
    gamma::Float64, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I,  A_matrix::AbstractMatrix{Float64}, 
      my_method, my_cache, row_vector::AbstractVector{Float64},vector_cache)
        # initialize the vectors
        num_lineages = 2
        row_vector .= 0
        vector_cache .= 0
        # find the index of the reverse_time before the est time 
        is_between = (reverse_comp_times .> start_time) .& (reverse_comp_times .<= coal_time)
        subset_times = reverse_comp_times[is_between]
        subset_E = reverse_E[is_between]
        subset_I = reverse_I[is_between]
        if length(subset_times) == 1
            # we're done
            E = subset_E[1]
            compI = subset_I[1]
            update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, compI)
            delta_t = coal_time - start_time
            L_vec = zeros(3)
            L_vec[2] = -(A_matrix[2,2] + A_matrix[2,3]) 
            exponential!(rmul!(A_matrix,delta_t), ExpMethodHigham2005(), my_cache)
            prob = transpose(A_matrix[start_state,:]) * L_vec
        elseif length(subset_times) > 1
            # then we need to do some things 
            comp_idx = 1
            E = subset_E[comp_idx]
            compI = subset_I[comp_idx]
            update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, compI)
            delta_t = subset_times[1] - start_time
            exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
            row_vector .= @view A_matrix[start_state,:]
            # this first event must have been alpha time, so update alpha_idx
            alpha_idx = findfirst(isequal(alpha_t), reverse_alpha_vec)
            alpha_idx += 1
            for t in 2:length(subset_times)
                alpha_t = reverse_alpha_vec[alpha_idx]
                E = subset_E[t]
                compI = subset_I[t]
                update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, compI)
                delta_t = subset_times[t] - subset_times[t-1]
                if t == length(subset_times)
                    # this is a coal_time
                    L_vec = zeros(3)
                    L_vec[2] = -(A_matrix[2,2] + A_matrix[2,3]) 
                    exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
                    #tranpose(expAt) * row_vector
                    mul!(vector_cache, transpose(A_matrix), row_vector)
                    # transpose(L) * row_vector
                    prob = transpose(vector_cache) * L_vec
                elseif t < length(subset_times)
                    exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
                    mul!(vector_cache, transpose(A_matrix), row_vector)
                    row_vector .= vector_cache
                    # this was another alpha event tick up alpha
                    alpha_idx += 1
                end 
            end 
        end 
        return prob
end 



"""
calc_samp_pdf_noalpha
"""
function calc_samp_pdf_noalpha(num_lineages::Int, end_time::Float64, start_state, start_time::Float64, 
    gamma::Float64, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I,  A_matrix::AbstractMatrix{Float64}, 
     row_vector::AbstractVector{Float64}, my_method, my_cache, my_output::AbstractVector{Float64}) 
    row_vector .= 0
    my_output .= 0
    is_between = (reverse_comp_times .> start_time) .& (reverse_comp_times .<= end_time)
    subset_times = reverse_comp_times[is_between]
    subset_E = reverse_E[is_between]
    subset_I = reverse_I[is_between]
    if length(subset_times) == 1
        E = subset_E[1]
        compI = subset_I[1]
        update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, compI)
        delta_t = end_time - start_time
        # A_matrix now holds the exponential of the A matrix
        exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
        my_output .= @view A_matrix[start_state, :]
    elseif length(subset_times) > 1
        # there is at least one alpha change time 
        comp_idx = 1
        E = subset_E[comp_idx]
        compI = subset_I[comp_idx]
        update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, compI)
        # update_L_matrix!(L_matrix, num_lineages, A_matrix)
        delta_t = subset_times[1] - start_time
        exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
        row_vector .= @view A_matrix[start_state,:]
        # this first event must have been alpha time, so update alpha_idx
        alpha_idx = findfirst(isequal(alpha_t), reverse_alpha_vec)
        alpha_idx += 1
        for t in 2:length(subset_times)
            alpha_t = reverse_alpha_vec[alpha_idx]
            E = subset_E[t]
            compI = subset_I[t]
            update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, compI)
            delta_t = subset_times[t] - subset_times[t-1]
            if t == length(subset_times)
                # this is a sample time
                exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
                #tranpose(expAt) * row_vector
                mul!(my_output, transpose(A_matrix), row_vector)
            elseif t < length(subset_times)
                exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
                mul!(my_output, transpose(A_matrix), row_vector)
                row_vector .= my_output
                # this was another alpha event tick up alpha
                alpha_idx += 1
            end 
        end
    end 
    return my_output, alpha_t
end  
"""
sample_internal_nodes_noalpha!
assumes that est_times includes alpha_times, coal_times and samp_times
the returned vector will not include the state at the time of last sampling, which is always known
"""
function sample_internal_nodes_noalpha!(num_lineages::Int, est_times::AbstractVector{Float64}, 
    coal_times::AbstractVector{Float64}, reverse_samp_times::AbstractVector, 
    reverse_samp_lin::AbstractVector, gamma::Float64, reverse_alpha_vec::AbstractVector{Float64}, 
    reverse_E::AbstractVector{Float64}, reverse_I::AbstractVector{Float64}, cache_dict::Dict{Int, Any}, mat_size::Int, 
    ks_dict::Dict{Int, Any},  expv_cache_dict::Dict{Int, Any}, A_matrix::AbstractMatrix{Float64}, 
    L_matrix::AbstractMatrix{Float64},my_output, vector_cache, row_vector, sampled_node_states::Vector{Int}, tstep_cutoff)
    # initialize many things
    # this is only for debugging
    # ll_vec = zeros(length(est_times))
    # initialize the active lineages and start and est times
    my_method = ExpMethodHigham2005()
    active_lineages = num_lineages
    active_start_time = 0.0
    active_est_time = est_times[1]
    # reset sampled_node_states
    sampled_node_states .= 0
    sampled_node_states[end] = 1
    start_state = 1
# set the ticker for sampling vs coalescent times vs alpha times 
    coal_time_idx = 1
    next_coal_time = coal_times[coal_time_idx]
    samp_time_idx = 1
    next_samp_time = if length(reverse_samp_times) > 0
        reverse_samp_times[samp_time_idx]
    else 
        -1
    end 
    # set the initial values of alpha, E and I
    alpha_t = reverse_alpha_vec[1]
    # the first value of reverse_E and reverse_I corresponds to the time before the time of last sampling 
    # so the first value corresponds to the first estimation time
    log_lik = 0.0
    # set initial distribution to 0 and 1 at the first state
    @views for t in 1:(length(est_times)-1)
        old_alpha_t = alpha_t
        delta_t = active_est_time - active_start_time
        if active_est_time == next_coal_time 
            # print("t")
            # print(t)
            # print("alpha_t")
            # print(alpha_t)
            active_pdf, alpha_t = calc_coal_pdf_noalpha(active_lineages, active_est_time, start_state, active_start_time, 
            gamma, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I,  
            A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
            L_matrix[1:(active_lineages + 1), 1:(active_lineages - 1)], 
            my_output[1:(active_lineages - 1)], my_method, cache_dict[(active_lineages + 1)], 
            row_vector[1:(active_lineages + 1)], vector_cache[1:(active_lineages + 1)]) 
            # sample proportional to the pdf 
            state = sample(1:(active_lineages-1), Weights(active_pdf))
            sampled_node_states[t] = state
            # we need to check if we're over the krylov size,
            # if we are, we need to recalculate the pdf using krylov, otherwise it won't match 
            if active_lineages + 1 > mat_size && delta_t < tstep_cutoff
                new_ll, alpha_t = calc_coal_ll_noalpha_krylov(A_matrix[1:(active_lineages+1),1:(active_lineages+1)], 
    L_matrix[1:(active_lineages+1), 1:(active_lineages-1)], my_vec[1:(active_lineages+1)],  ks_dict[active_lineages + 1], expv_cache_dict[active_lineages + 1], 
    active_lineages, active_est_time, state, active_start_time, start_state,
    gamma, old_alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I)
                active_start_state = sampled_node_states[t]
                # CHECK THESE
                log_lik += new_ll
            else 
                # we can just use what we have already calculated
                log_lik += log(active_pdf[state])
                # ll_vec[t] = log(active_pdf[state])
            end
            coal_time_idx += 1
            next_coal_time = coal_times[coal_time_idx]
            active_lineages -= 1
        elseif active_est_time == next_samp_time 
            pmf, alpha_t = calc_samp_pdf_noalpha(active_lineages, active_est_time, start_state, active_start_time, 
            gamma, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I,  A_matrix[1:(active_lineages+1), 1:(active_lineages+1)], 
            row_vector[1:(active_lineages+1)], my_method, cache_dict[(active_lineages+1)], my_output[1:(active_lineages+1)]) 
            state = sample(1:(active_lineages + 1), Weights(pmf))
            sampled_node_states[t] = state
            if active_lineages + 1 > mat_size && delta_t < tstep_cutoff
                new_ll, alpha_t = calc_samp_ll_noalpha_krylov(A_matrix[1:(active_lineages+1), 1:(active_lineages+1)], 
     ks_dict[(active_lineages+1)], expv_cache_dict[(active_lineages+1)], 
    active_lineages, active_est_time, state, active_start_time, start_state,
    gamma, old_alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I)
                log_lik += new_ll
            else 
                log_lik += log(pmf[state])
            end
            active_lineages += reverse_samp_lin[samp_time_idx]
            samp_time_idx += 1
            if samp_time_idx <= length(reverse_samp_times)
                next_samp_time = reverse_samp_times[samp_time_idx]
            end
        end
        # update the active init dist
        start_state = state
        # update the active start time
        active_start_time = active_est_time
        # update the active time
        active_est_time = est_times[t + 1]
    end
    last_prob = calc_coal_twolin_pdf_noalpha(active_est_time, start_state, active_start_time, 
    gamma::Float64, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I,  A_matrix[1:3, 1:3], 
     my_method, cache_dict[3], row_vector[1:3],
     vector_cache[1:3])
    log_lik += log(last_prob)
    # ll_vec[end] = log(last_prob)
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
    return log_lik, sampled_node_states
end 
