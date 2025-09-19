using DataFrames
"""
Julia version of ei_inv_event_rate.R 
"""
function ei_event_rate_sampler(t::Float64; 
                            init_eff_frame::DataFrame, 
                            alpha::Float64, 
                            gamma::Float64, 
                            max_time::Float64, 
                            nE::Int64, 
                            nI::Int64)

    if t > max_time
        rate = Inf
        coal_rate = Inf
        I2E_rate = 0
        E2I_rate = 0
        E = 0
        I = 1
        return rate, coal_rate, I2E_rate, E2I_rate, E, I
    else
        # Find the index of reverse_time closest to t[i]
        rate_idx = findlast(x -> x >= t, init_eff_frame.reverse_time)

        I = init_eff_frame.I[rate_idx]
        E = init_eff_frame.E[rate_idx]

        # Calculate the rates with legal move indicators
        coal_rate = (alpha * nE * nI / (E + 1E-10)) 
        I2E_rate = (gamma * nI * (E + 1) / (I + 1E-10)) 
        E2I_rate = max(0, (alpha * nE * (I - nI) / (E + 1E-10))) 

        rate = coal_rate + I2E_rate + E2I_rate
    end
    return rate, coal_rate, I2E_rate, E2I_rate, E, I
end