"""
calc_coal_ll_noalpha_higham()
For known start and times and states, calculate the likelihood contribution of a coalescent interval
using full matrix exponentiation using scaling and squaring from Higham (2005)

returns log likelihood contribution and value of alpha at end_time
# Arguments 
-num_lineages: number of extant lineages
-coal_time: time of coalescence 
-end_state: ending coalescent state
-start_time: interval start time
-start_state : start state
-gamma: gamma param
-alpha_t : value of alpha at start_time
-reverse_alpha_vec : vector of alpha values in order of reverse_time
-reverse_comp_times: vector of reverse times in which E and I change values
-reverse_E: vector of values of E at times reverse_comp_times (same length as reverse_comp_times)
-reverse_E: vector of values of I at times reverse_comp_times (same length as reverse_comp_times)
-A_matrix: matrix of rates of transitioning to non-coalescent states (n+1 X n+1)
-L_matrix: matrix of rates of transitioning to coalescent states (n+ 1 X n-1)
-my_method: method for matrix exponentiation, set to ExpMethodHigham2005()
-my_cache: cache for matrix exponentiation
-row_vector: place to store matrix multiplication 
-vector_cache: other place to store matrix multiplication
"""
function calc_coal_ll_noalpha_higham(num_lineages::Int, coal_time::Float64, end_state::Int64, start_time::Float64,  
    start_state::Int64, gamma::Float64, alpha_t::Float64, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I,  
    A_matrix::AbstractMatrix{Float64}, L_matrix::AbstractMatrix{Float64}, my_method, my_cache, 
    row_vector::AbstractVector{Float64}, vector_cache) 
    # initialize the vectors
    row_vector .= 0
    vector_cache .= 0
    # find the index of the reverse_time before the est time 
    is_between = (reverse_comp_times .> start_time) .& (reverse_comp_times .<= coal_time)
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
        delta_t = coal_time - start_time
        # A_matrix now holds the exponential of the A matrix
        exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
        # this is the product term (exp(At)
        ll_prob = transpose(L_matrix[:, end_state]) * A_matrix[start_state, :]
    elseif length(subset_times) > 1
        # now we need to do something more complicated
        # there is at least one alpha update in here
        # go in forward time b/c it means less matrix matrix
        comp_idx = 1
        E = subset_E[comp_idx]
        I = subset_I[comp_idx]
        update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, I)
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
            I = subset_I[t]
            update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, I)
            delta_t = subset_times[t] - subset_times[t-1]
            if t == length(subset_times)
                # this is a coal_time
                update_L_matrix!(L_matrix, num_lineages, A_matrix)
                exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
                #tranpose(expAt) * row_vector
                mul!(vector_cache, transpose(A_matrix), row_vector)
                # transpose(L) * row_vector
                ll_prob = transpose(L_matrix[:, end_state]) * vector_cache
            elseif t < length(subset_times)
                exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
                mul!(vector_cache, transpose(A_matrix), row_vector)
                row_vector .= vector_cache
                # this was another alpha event tick up alpha
                alpha_idx += 1
            end 
        end 
    end 
    return log(ll_prob), alpha_t
end 
"""
calc_samp_ll_noalpha_higham()
For known start and times and states, calculate the likelihood contribution of a sampling interval
using full matrix exponentiation using scaling and squaring from Higham (2005)

returns log likelihood contribution and value of alpha at end_time
# Arguments 
-num_lineages: number of extant lineages
-end_state: ending coalescent state
-end_time: time of the sampling event
-start_state : start state
-start_time: interval start time
-gamma: gamma param
-alpha_t : value of alpha at start_time
-reverse_alpha_vec : vector of alpha values in order of reverse_time
-reverse_comp_times: vector of reverse times in which E and I change values
-reverse_E: vector of values of E at times reverse_comp_times (same length as reverse_comp_times)
-reverse_E: vector of values of I at times reverse_comp_times (same length as reverse_comp_times)
-A_matrix: matrix of rates of transitioning to non-coalescent states (n+1 X n+1)
-row_vector: place to store matrix multiplication 
-my_method: method for matrix exponentiation, set to ExpMethodHigham2005()
-my_cache: cache for matrix exponentiation
-vector_cache: other place to store matrix multiplication
"""
function calc_samp_ll_noalpha_higham(num_lineages::Int,  end_state, end_time::Float64, start_state, start_time::Float64, 
    gamma::Float64, alpha_t, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I,  A_matrix::AbstractMatrix{Float64}, 
     row_vector::AbstractVector{Float64}, my_method, my_cache, vector_cache::AbstractVector{Float64}) 
    row_vector .= 0
    vector_cache .= 0
    is_between = (reverse_comp_times .> start_time) .& (reverse_comp_times .<= end_time)
    subset_times = reverse_comp_times[is_between]
    subset_E = reverse_E[is_between]
    subset_I = reverse_I[is_between]
    if length(subset_times) == 1
        E = subset_E[1]
        I = subset_I[1]
        update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, I)
        delta_t = end_time - start_time
        # A_matrix now holds the exponential of the A matrix
        exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
        ll_prob = A_matrix[start_state, end_state]
    elseif length(subset_times) > 1
        # there is at least one alpha change time 
        comp_idx = 1
        E = subset_E[comp_idx]
        I = subset_I[comp_idx]
        update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, I)
        # update_L_matrix!(L_matrix, num_lineages, A_matrix)
        delta_t = subset_times[1] - start_time
        exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
        row_vector .= @view  A_matrix[start_state,:]
        # this first event must have been alpha time, so update alpha_idx
        alpha_idx = findfirst(isequal(alpha_t), reverse_alpha_vec)
        alpha_idx += 1
        for t in 2:length(subset_times)
            alpha_t = reverse_alpha_vec[alpha_idx]
            E = subset_E[t]
            I = subset_I[t]
            update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, I)
            delta_t = subset_times[t] - subset_times[t-1]
            if t == length(subset_times)
                # this is a sample time
                exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
                #tranpose(expAt) * row_vector
                # CHECK THIS
                ll_prob =  transpose(A_matrix[:, end_state]) * row_vector
            elseif t < length(subset_times)
                exponential!(rmul!(A_matrix,delta_t), my_method, my_cache)
                mul!(vector_cache, transpose(A_matrix), row_vector)
                row_vector .= vector_cache
                # this was another alpha event tick up alpha
                alpha_idx += 1
            end 
        end
    end 
    return log(ll_prob), alpha_t
end  
