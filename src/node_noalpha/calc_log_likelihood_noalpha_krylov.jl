"""
for known start and end times and states
    calculate pdf of coalescent interval with possible alpha changes
    using krylov subpsace method
"""
function calc_coal_ll_noalpha_krylov(A_matrix::AbstractMatrix{Float64}, 
    L_matrix::AbstractMatrix{Float64}, my_vec::AbstractVector{Float64}, my_ks, expv_cache, 
    active_lineages::Int, end_time::Float64, end_state::Int, start_time::Float64, start_state,
    gamma::Float64, alpha_t::Float64, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I)

    my_output .= 0
    row_vector .= 0
    vector_cache .= 0
    # find the index of the reverse_time before the est time 
    is_between = (reverse_comp_times .> start_time) .& (reverse_comp_times .<= end_time)
    subset_times = reverse_comp_times[is_between]
    subset_E = reverse_E[is_between]
    subset_I = reverse_I[is_between]
    if length(subset_times) == 1
        E = subset_E[1]
        I = subset_I[1]
        update_A_matrix_cm_simplev2!(A_matrix, active_lineages, alpha_t, gamma, E, I)
        update_L_matrix!(L_matrix, active_lineages, A_matrix)
        delta_t = end_time - start_time
        my_vec .= @view L_matrix[:,end_state]
        ExponentialUtilities.arnoldi!(my_ks, A_matrix, my_vec; ishermitian=false)
        ExponentialUtilities.expv!(my_vec, delta_t, my_ks, cache = expv_cache)
        ll_prob = my_vec[start_state]
    elseif length(subset_times) > 1
        # more stuff happens, there is at least one alpha change 
        comp_idx = 1
        E = subset_E[comp_idx]
        I = subset_I[comp_idx]
        index_vector = zeros(active_lineages + 1)
        index_vector[start_state] = 1.0
        update_A_matrix_cm_simplev2!(A_matrix, active_lineages, alpha_t, gamma, E, I)
        delta_t = subset_times[1] - start_time
        # tomorrow we need to check that this works
        ExponentialUtilities.arnoldi!(my_ks, transpose(A_matrix), index_vector; ishermitian=false)
        ExponentialUtilities.expv!(index_vector, delta_t, my_ks, cache = expv_cache)
        # this first event must have been alpha time, so update alpha_idx
        alpha_idx = findfirst(isequal(alpha_t), reverse_alpha_vec)
        alpha_idx += 1
        for t in 2:length(subset_times)
            alpha_t = reverse_alpha_vec[alpha_idx]
            E = subset_E[t]
            I = subset_I[t]
            update_A_matrix_cm_simplev2!(A_matrix, active_lineages, alpha_t, gamma, E, I)
            delta_t = subset_times[t] - subset_times[t-1]
            if t == length(subset_times)
                # this is a coal_time
                update_L_matrix!(L_matrix, active_lineages, A_matrix)
                ExponentialUtilities.arnoldi!(my_ks, transpose(A_matrix), index_vector; ishermitian=false)
                ExponentialUtilities.expv!(index_vector, delta_t, my_ks, cache = expv_cache)
                #tranpose(expAt) * row_vector
                # transpose(L_vector) * row_vector
                ll_prob = transpose(L_matrix[:, end_state]) * index_vector
            elseif t < length(subset_times)
                ExponentialUtilities.arnoldi!(my_ks, transpose(A_matrix), index_vector; ishermitian=false)
                ExponentialUtilities.expv!(index_vector, delta_t, my_ks, cache = expv_cache)
                # this was another alpha event tick up alpha
                alpha_idx += 1
            end 
        end
    end
    return log(ll_prob), alpha_t
end   
"""
calc_samp_pdf_noalpha_krylov
calculate ll contribution of sampling event with known start and end states
and assuming any alpha changes get rolled in
"""
function calc_samp_ll_noalpha_krylov(A_matrix::AbstractMatrix{Float64}, 
     my_ks, expv_cache, 
    active_lineages::Int, end_time::Float64, end_state::Int, start_time::Float64, start_state,
    gamma::Float64, alpha_t::Float64, reverse_alpha_vec, reverse_comp_times, reverse_E, reverse_I)
    my_output .= 0
    row_vector .= 0
    vector_cache .= 0
    # find the index of the reverse_time before the est time 
    is_between = (reverse_comp_times .> start_time) .& (reverse_comp_times .<= end_time)
    subset_times = reverse_comp_times[is_between]
    subset_E = reverse_E[is_between]
    subset_I = reverse_I[is_between]
    if length(subset_times) == 1
        E = subset_E[1]
        I = subset_I[1]
        update_A_matrix_cm_simplev2!(A_matrix, active_lineages, alpha_t, gamma, E, I)
        delta_t = end_time - start_time
        my_vec = zeros(active_lineages + 1)
        my_vec[start_state] = 1.0
        ExponentialUtilities.arnoldi!(my_ks, transpose(A_matrix), my_vec; ishermitian=false)
        ExponentialUtilities.expv!(my_vec, delta_t, my_ks, cache = expv_cache)
        ll_prob = my_vec[end_state]
    elseif length(subset_times) > 1
        # more stuff happens, there is at least one alpha change 
        comp_idx = 1
        E = subset_E[comp_idx]
        I = subset_I[comp_idx]
        index_vector = zeros(active_lineages + 1)
        index_vector[start_state] = 1.0
        update_A_matrix_cm_simplev2!(A_matrix, active_lineages, alpha_t, gamma, E, I)
        delta_t = subset_times[1] - start_time
        # tomorrow we need to check that this works
        ExponentialUtilities.arnoldi!(my_ks, transpose(A_matrix), index_vector; ishermitian=false)
        ExponentialUtilities.expv!(index_vector, delta_t, my_ks, cache = expv_cache)
        # this first event must have been alpha time, so update alpha_idx
        alpha_idx = findfirst(isequal(alpha_t), reverse_alpha_vec)
        alpha_idx += 1
        for t in 2:length(subset_times)
            alpha_t = reverse_alpha_vec[alpha_idx]
            E = subset_E[t]
            I = subset_I[t]
            update_A_matrix_cm_simplev2!(A_matrix, active_lineages, alpha_t, gamma, E, I)
            delta_t = subset_times[t] - subset_times[t-1]
            if t == length(subset_times)
                # this is a samp time and we're done
                ExponentialUtilities.arnoldi!(my_ks, transpose(A_matrix), index_vector; ishermitian=false)
                ExponentialUtilities.expv!(index_vector, delta_t, my_ks, cache = expv_cache)
                #tranpose(expAt) * row_vector
                # transpose(L_vector) * row_vector
                ll_prob = index_vector[end_state]
            elseif t < length(subset_times)
                ExponentialUtilities.arnoldi!(my_ks, transpose(A_matrix), index_vector; ishermitian=false)
                ExponentialUtilities.expv!(index_vector, delta_t, my_ks, cache = expv_cache)
                # this was another alpha event tick up alpha
                alpha_idx += 1
            end 
        end
    end
    return log(ll_prob), alpha_t
end   


