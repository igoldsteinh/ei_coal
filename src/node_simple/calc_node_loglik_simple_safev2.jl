# this file is a rewrite of calc_node_loglik_joint_hetchron.jl
# with the added assumption that the rates are constant
# between estimation states
# this removes the need for much indexing and for loops
"""
calc_node_loglik_simple_safe
should calculate node log likelihoods under the assumption that 
the rates only change at the estimation times
not too much trickery with re-use of variables
there will be a revised version for use in the model 
"""
function calc_node_loglik_simple_safev2(num_lineages, est_times, coal_times, est_states, start_time, last_samp_time, 
    reverse_samp_times, reverse_samp_lin, alpha_times, gamma, alpha_vec, reverse_E, reverse_I, total_pop,
 mat_size, lin_E, lin_I)
    # construct needed things 
    ll_vec = Vector{Float64}(undef, length(est_times))
    lik_vec = Vector{Vector{Float64}}(undef, length(est_times))
    sizes = 2:(max_lineages + 1)
    my_method = ExpMethodHigham2005()
    cache_dict = preallocate_caches(sizes, my_method)
    ks_dict = preallocate_krylov((mat_size+1):max_lineages + 1)
    expv_cache_dict = preallocate_expv_cache((mat_size+1):max_lineages + 1)
    A_matrix = zeros( max_lineages+1, max_lineages+1)
    L_matrix = zeros( max_lineages+1, max_lineages-1)
    my_method = ExpMethodHigham2005()
    L_vec = zeros(3)
    my_vec = zeros(max_lineages + 1)
    temp = similar(my_vec)
    pop_big_enough = zeros(length(curr_lin))



    # likelihood checks, pop size too large, or if I is too small
    if any(total_pop .> 8E9)
        return -Inf, ll_vec, lik_vec
    end 
    # check if the populations are big enough at each coal and sample time
    check_total_pop_sizev2!(lin_E, lin_I, reverse_E, reverse_I, pop_big_enough)
    # if they're ever not, reject immediately
    if any(x -> x == 0, pop_big_enough)
        return -Inf, ll_vec, lik_vec
    end
    # set up active_lineages, active_states, and active_times
    # these are the starting times and states, and the ending times and states for each likelihood chunk
    # they are updated as we move through estimation times
    active_lineages = num_lineages
    active_start_time = start_time
    active_est_time = est_times[1]
    est_state_idx = 1
    active_start_state = est_states[est_state_idx]
    active_est_state = est_states[est_state_idx + 1]
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
    ll_vec = Vector{Float64}(undef, length(est_times))
    lik_vec = Vector{Vector{Float64}}(undef, length(est_times))
    delta_t = active_est_time - active_start_time
    @inbounds for j in 1:(length(est_times) - 1)
        if active_est_time == next_coal_time
            if (active_lineages + 1) > mat_size && delta_t < 1.0 # use the krylov method
                @views lik_vec[j] = calc_coal_pdf_vecform_krylov_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                L_matrix[1:(active_lineages + 1), 1:(active_lineages - 1)], my_vec[1:(active_lineages + 1)], 
                ks_dict[Int(active_lineages + 1)], expv_cache_dict[Int(active_lineages + 1)], active_lineages, active_est_time, 
                active_est_state, active_start_time, gamma, alpha_t, E, I)
            else # use the matrix method
                @views lik_vec[j] = calc_coal_pdf_vecform_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                L_matrix[1:(active_lineages + 1), 1:(active_lineages - 1)], my_vec[1:(active_lineages + 1)], my_method, 
                cache_dict[active_lineages + 1], active_lineages, active_est_time, active_est_state, 
                active_start_time, gamma, alpha_t, E, I)
            end
            ll_vec[j] = if lik_vec[j][active_start_state] > 0
                log(lik_vec[j][active_start_state])
            else
                -Inf
            end
            log_lik += if lik_vec[j][active_start_state] > 0
                log(lik_vec[j][active_start_state])
            else
                -Inf
            end
            active_lineages -= 1
            coal_time_idx += 1
            next_coal_time = coal_times[coal_time_idx]
        elseif active_est_time == next_samp_time
            if (active_lineages + 1) > mat_size && delta_t < 1.0 # use the krylov method 
                @views lik_vec[j] = calc_samp_lik_vecform_krylov_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                my_vec[1:(active_lineages + 1)], ks_dict[Int(active_lineages + 1)], expv_cache_dict[Int(active_lineages + 1)], active_lineages, 
                active_est_time, active_est_state, active_start_time, gamma, alpha_t, E, I)
            else # use the matrix method
                @views lik_vec[j] = calc_samp_lik_vecform_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                my_vec[1:(active_lineages + 1)], my_method, cache_dict[active_lineages + 1], temp[1:(active_lineages + 1)], 
                active_lineages, active_est_time, active_est_state, active_start_time, gamma, alpha_t, E, I)
            end 
            ll_vec[j] = log(lik_vec[j][active_start_state])
            log_lik += log(lik_vec[j][active_start_state])    
            active_lineages += reverse_samp_lin[samp_time_idx]
            samp_time_idx += 1
            if samp_time_idx <= length(reverse_samp_times)
                next_samp_time = reverse_samp_times[samp_time_idx]
            end
        elseif active_est_time == next_alpha_time
            # calculate alpha time using sampling functions
            # do not update active lineages 
            # but do update alpha_t (after likleihood calculation)
            if (active_lineages + 1) > mat_size && delta_t < 1.0 # use the krylov method 
                @views lik_vec[j] = calc_samp_lik_vecform_krylov_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                my_vec[1:(active_lineages + 1)], ks_dict[Int(active_lineages + 1)], expv_cache_dict[Int(active_lineages + 1)], active_lineages, 
                active_est_time, active_est_state, active_start_time, gamma, alpha_t, E, I)
            else # use the matrix method
                @views lik_vec[j] = calc_samp_lik_vecform_simplev2(A_matrix[1:(active_lineages + 1), 1:(active_lineages + 1)], 
                my_vec[1:(active_lineages + 1)], my_method, cache_dict[active_lineages + 1], temp[1:(active_lineages + 1)], 
                active_lineages, active_est_time, active_est_state, active_start_time, gamma, alpha_t, E, I)
            end 
            ll_vec[j] = log(lik_vec[j][active_start_state])
            log_lik += log(lik_vec[j][active_start_state])    
            # update alpha_t
            alpha_time_idx -= 1
            if alpha_time_idx > 0
                next_alpha_time = abs(alpha_times[alpha_time_idx] - last_samp_time)
                alpha_t = alpha_vec[alpha_time_idx]
            else
                next_alpha_time = -1
            end 
        end
        # update the start and end states
        est_state_idx += 1
        active_start_state = est_states[est_state_idx]
        active_est_state = est_states[est_state_idx + 1]
        # update the active start time
        active_start_time = active_est_time
        # update the active est time
        active_est_time = est_times[j + 1]
        # update E and I 
        comp_idx += 1
        E = reverse_E[comp_idx]
        I = reverse_I[comp_idx]
        delta_t = active_est_time - active_start_time
    end 
    # now do the the two lineage
    @views lik_vec[end] = calc_twolin_pdf_vecform_simplev2(A_matrix[1:3, 1:3], L_vec, 
    my_vec[1:3], my_method, cache_dict[3], active_est_time, active_start_time, gamma, alpha_t, E, I)
    ll_vec[end] = log(lik_vec[end][active_start_state])
    log_lik += log(lik_vec[end][active_start_state])
    return log_lik, ll_vec, lik_vec
end 

"""
calc_samp_lik_vecform_krylov_simple
this is the likelihood contribution of a sampling event in vector form, 
choose the correct entry of the column based on the start state 
note this should work for sampling events and also for alpha changes 
calculated with krylov subspace methods
"""
function calc_samp_lik_vecform_krylov_simplev2(A_matrix, my_vec, my_ks, expv_cache, active_lineages, est_time, branch_state, start_time, 
    gamma, alpha_t, E, I)
    my_vec .= 0
    my_vec[branch_state] = 1.0
    update_A_matrix_cm_simplev2!(A_matrix, active_lineages, alpha_t, gamma, E, I)
    delta_t = est_time - start_time
    # calculate action of exponential vector
    ExponentialUtilities.arnoldi!(my_ks, A_matrix, my_vec; ishermitian=false)
    ExponentialUtilities.expv!(my_vec, delta_t, my_ks, cache = expv_cache)
    return my_vec
end 

"""
calc_samp_lik_vecform_simple
this is the likelihood contribution of a sampling event in vector form, 
choose the correct entry of the column based on the start state 
note this should work for sampling events and also for alpha changes 
"""
function calc_samp_lik_vecform_simplev2(A_matrix, my_vec, my_method, my_cache, temp, active_lineages, est_time, branch_state, 
    start_time, gamma, alpha_t, E, I)
    my_vec .= 0.0
    my_vec[branch_state] = 1.0
    update_A_matrix_cm_simplev2!(A_matrix, active_lineages, alpha_t, gamma, E, I)
    delta_t = est_time - start_time
    mul!(temp, exponential!(rmul!(A_matrix, delta_t), my_method, my_cache), my_vec)
    my_vec .= temp
    return my_vec
end 
"""
calc_coal_pdf_vecform_simple
"""
function calc_coal_pdf_vecform_simplev2(A_matrix, L_matrix, my_vec, my_method, my_cache, num_lineages, coal_time, 
    coal_type, start_time, gamma, alpha_t, E, I)
    update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, I)
    update_L_matrix!(L_matrix, num_lineages, A_matrix)
    delta_t = coal_time - start_time
    mul!(my_vec, exponential!(rmul!(A_matrix,delta_t), my_method, my_cache), @view L_matrix[:,coal_type])
    return my_vec
end 
"""
calc_coal_pdf_vecform_krylov_simple
but we're trying to just do vector multiplication
"""
function calc_coal_pdf_vecform_krylov_simplev2(A_matrix, L_matrix, my_vec, my_ks, expv_cache, num_lineages, coal_time, coal_type, start_time, 
    gamma, alpha_t, E, I)
    update_A_matrix_cm_simplev2!(A_matrix, num_lineages, alpha_t, gamma, E, I)
    update_L_matrix!(L_matrix, num_lineages, A_matrix)
    # calculate first action of exponential vector
    delta_t = coal_time - start_time
    my_vec .= @view L_matrix[:,coal_type]
    if any(isnan.(A_matrix))|| any(isinf.(A_matrix)) || any(isnan.(my_vec)) || any(isinf.(my_vec))
        my_vec .= 0.0
        return my_vec
    end
    ExponentialUtilities.arnoldi!(my_ks, A_matrix, my_vec; ishermitian=false)
    ExponentialUtilities.expv!(my_vec, delta_t, my_ks, cache = expv_cache)
    return my_vec
end 
"""
calc_twolin_pdf_vecform_simple
caculate the pdf of a two lineage coalescent even with known stard and end states
vecform version of two lineage pdf
WHY ARE YOU ALLOCATING SO MUCH?
"""
function calc_twolin_pdf_vecform_simplev2(A_matrix, L_vec, my_vec, my_method, my_cache,  coal_time,  start_time, 
    gamma, alpha_t, E, I)
        # find the index of the reverse_time before the coalescent time 
        active_lineages = 2
        update_A_matrix_cm_simplev2!(A_matrix, active_lineages, alpha_t, gamma, E, I)
        delta_t = coal_time - start_time
        L_vec.=0.0
        L_vec[2] = -(A_matrix[2,2] + A_matrix[2,3]) 
        mul!(my_vec, exponential!(rmul!(A_matrix,delta_t), my_method, my_cache), L_vec)
        return my_vec
end 
"""
column major update_A_matrix!
this seems to be slighly slower than the other one
hopefully the column major has better performance in later functions
but we could experiment if needed 
num_lineages is the current number of extant lineages
E, I my_alpha, gamma are the current values of alpha, gamma, E and I for the rates 
we are going to do something odd and remove all indicators enforcing lineage size against compartment size
we will make it so that the E2I rate is 0 if nI < I 
"""
function update_A_matrix_cm_simplev2!(A_matrix, num_lineages, my_alpha, gamma, E, I)
    # first set all entries to 0
    A_matrix .= 0.0
    # Precompute constants to avoid redundant calculations
    # Iterate over columns instead of rows
    # column 1
    I_lineages = num_lineages
    E_lineages = num_lineages - I_lineages
    # precompute constants?
    inv_I = 1.0 / (I + 1E-10)
    inv_E = 1.0 / (E + 1E-10)
    
    # A_matrix[1,1] = -(I_lineages * gamma * (E + 1) / (I + 1E-10) * (I_lineages <= I ? 1.0 : 0.0) * (E_lineages <= E ? 1.0 : 0.0))
    I2E_rate = I_lineages * gamma * (E + 1) * inv_I 
    A_matrix[1,1] = -I2E_rate
    # one down is E2I from I_lineages - 1 E_lineages + 1
    E2I_rate = max(0.0, (E_lineages + 1) * my_alpha * (I - (I_lineages - 1)) * inv_E)
    A_matrix[2,1] = E2I_rate
    # A_matrix[2,1] = (E_lineages + 1) * my_alpha * (I - (I_lineages - 1)) * inv_E 
     
    @inbounds for j in 2:(num_lineages)
            I_lineages = num_lineages - j + 1
            E_lineages = num_lineages - I_lineages
            # now update one above the diagonal, this is an I2E move from I_lineages + 1 E_lineages - 1
            # A_matrix[j - 1,j] = (I_lineages + 1) * gamma * (E + 1) * inv_I * ((I_lineages + 1) <= I ? 1.0 : 0.0) * ((E_lineages - 1) <= E ? 1.0 : 0.0)
            I2E_rate = (I_lineages + 1) * gamma * (E + 1) * inv_I 
            A_matrix[j - 1,j] = I2E_rate

            # update the diagonal
            # A_matrix[j,j] = -(E_lineages * my_alpha * (I - I_lineages) * inv_E * (I_lineages + 1 <= I ? 1.0 : 0.0) * (E_lineages <= E ? 1.0 : 0.0) + 
            # I_lineages * gamma * (E + 1) * inv_I * (I_lineages <= I ? 1.0 : 0.0) * (E_lineages <= E ? 1.0 : 0.0) + 
            # (I_lineages * E_lineages * my_alpha * inv_E * (I_lineages <= I ? 1.0 : 0.0) * (E_lineages <= E ? 1.0 : 0.0)))
            total_rate = max(0, (E_lineages * my_alpha * (I - I_lineages) * inv_E) + 
            I_lineages * gamma * (E + 1) * inv_I  + 
            (I_lineages * E_lineages * my_alpha * inv_E))
            A_matrix[j,j] = -total_rate
            # now update one below the diagonal, this is an E2I move from I_lineages - 1 E_lineages + 1
            E2I_rate = max(0, (E_lineages + 1) * my_alpha * (I - (I_lineages - 1)) * inv_E)
            A_matrix[j + 1,j] = E2I_rate
            # A_matrix[j + 1,j] = (E_lineages + 1) * my_alpha * (I - (I_lineages - 1)) * inv_E 
    end
    # now update the last column 
    I_lineages = 0
    E_lineages = num_lineages
    j = num_lineages + 1
    # one above the diagonal is an I2E move 
    # A_matrix[j - 1, j] = (I_lineages + 1) * gamma * (E + 1) * inv_I * ((I_lineages + 1) <= I ? 1.0 : 0.0) * ((E_lineages - 1) <= E ? 1.0 : 0.0)
    I2E_rate = (I_lineages + 1) * gamma * (E + 1) * inv_I 
    A_matrix[j - 1, j] = I2E_rate
    # now do the diagonal
    # A_matrix[j,j] = -(E_lineages * my_alpha * (I - I_lineages) * inv_E * (I_lineages + 1 <= I ? 1.0 : 0.0) * (E_lineages <= E ? 1.0 : 0.0))
    E2I_rate = max(0, (E_lineages * my_alpha * (I - I_lineages) * inv_E))
    A_matrix[j,j] = -E2I_rate
end
"""
create_L_matrix
Creates matrix of rates to transition into absorbing states
with n lineages, there are n-1 absorbing states
"""
function update_L_matrix!(L_matrix, num_lineages, A_matrix)
    L_matrix .= 0.0
    @views for j in 1:(num_lineages -1)
        L_matrix[j + 1, j] = -sum(A_matrix[j + 1, :])
    end 
    return L_matrix
end 
"""
preallocate caches for all sizes of A matrices
"""
function preallocate_caches(sizes, method)
    caches = Dict{Int, Any}()
    for (cols) in sizes
        A_matrix = zeros(cols, cols)  # Temporary matrix for cache allocation
        caches[cols] = ExponentialUtilities.alloc_mem(A_matrix, method)
    end
    return caches
end
"""
preallocate krylove spaces for all sizes of A matrices
"""
function preallocate_krylov(sizes)
    k_space = Dict{Int, Any}()
    for (cols) in sizes
        state = zeros(cols)  # Temporary matrix for cache allocation
        krylovdim = 30 # default
        k_space[cols] = ExponentialUtilities.KrylovSubspace{Float64}(length(state), min(krylovdim, length(state)))
    end
    return k_space
end
"""
preallocate expv cache for all sizes of A matrices
"""
function preallocate_expv_cache(sizes)
    expv_space = Dict{Int, Any}()
    for (cols) in sizes
        expv_space[cols] = ExponentialUtilities.ExpvCache{Float64}(cols)
    end
    return expv_space
end



"""
check_comp_pop_size(lin_E, lin_I, reverse_E, reverse_I)
assuming reverse_E, reverse_I are same length as lin_E, lin_I
check that the populations in the compartments are large enough 
"""
function check_total_pop_sizev2!(lin_E, lin_I, reverse_E, reverse_I, pop_big_enough)
    # check that the total population size is large enough
    for i in 1:length(lin_E)
        if floor(reverse_E[i]) + floor(reverse_I[i]) < lin_E[i] + lin_I[i]
            pop_big_enough[i] = false
        else
            pop_big_enough[i] = true
        end
    end
    return pop_big_enough 
end
