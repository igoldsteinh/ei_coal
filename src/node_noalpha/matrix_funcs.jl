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
function update_A_matrix_cm_simplev2!(A_matrix::AbstractMatrix{Float64}, num_lineages::Int, my_alpha::Float64, gamma::Float64, E::Float64, I::Float64)
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
    E2I_rate = max(0, (E_lineages + 1) * my_alpha * (I - (I_lineages - 1)) * inv_E)
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
function update_L_matrix!(L_matrix::AbstractMatrix{Float64}, num_lineages::Int, A_matrix::AbstractMatrix{Float64})
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
