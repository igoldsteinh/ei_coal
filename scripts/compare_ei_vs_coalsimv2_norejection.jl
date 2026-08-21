# this is a file which simulates EI trajectories
# without any rejection condition on the final population size
# used to see the effect of the rejection condition (I[end] >= 5/6)
# on the distribution of stochastic trajectories, across R0 settings
using DrWatson
using CSV
using Distributions
using Statistics
using StatsBase
using Random
include(srcdir("sim_ei.jl"))

nu = 1/3
gamma = 1/2
stop_time = 35.0
I_init = 1
total_pop_sims = 100

### r0 = 2 ------------------------------------------------------------------
r0 = 2
test_alpha = r0 * nu
Random.seed!(235)
all_state_dfs = []
for i in 1:total_pop_sims
    _, state_frame = sim_ei(I_init, test_alpha, nu, gamma, stop_time)
    state_frame.sim .= i
    push!(all_state_dfs, state_frame)
end
final_state_frame = vcat(all_state_dfs...)
CSV.write(datadir("compare_data", "compare_ei_state_frame_norejection.csv"), final_state_frame)

### r0 = 1.5 ------------------------------------------------------------------
r0 = 1.5
test_alpha = r0 * nu
Random.seed!(2357)
all_state_dfs = []
for i in 1:total_pop_sims
    _, state_frame = sim_ei(I_init, test_alpha, nu, gamma, stop_time)
    state_frame.sim .= i
    push!(all_state_dfs, state_frame)
end
final_state_frame = vcat(all_state_dfs...)
CSV.write(datadir("compare_data", "compare_ei_state_frame_norejection_r=1.5.csv"), final_state_frame)

### r0 = 1.3 ------------------------------------------------------------------
r0 = 1.3
test_alpha = r0 * nu
Random.seed!(2357)
all_state_dfs = []
for i in 1:total_pop_sims
    _, state_frame = sim_ei(I_init, test_alpha, nu, gamma, stop_time)
    state_frame.sim .= i
    push!(all_state_dfs, state_frame)
end
final_state_frame = vcat(all_state_dfs...)
CSV.write(datadir("compare_data", "compare_ei_state_frame_norejection_r=1.3.csv"), final_state_frame)
