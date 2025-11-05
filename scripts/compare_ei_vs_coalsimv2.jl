# this is a file which simulates an EI trajectory
# constructs true EI coalescent trees
# and also simulates trees using simcoal
using DrWatson
using CSV
using Distributions
using Statistics
using StatsBase
using Combinatorics
using Random
using ExponentialUtilities
include(srcdir("sim_ei.jl"))
include(srcdir("create_coal_trees.jl"))
include(srcdir("time_transform", "ei_event_rate_sampler.jl"))
include(srcdir("time_transform", "sim_next_eicoal_time_tt.jl"))
include(srcdir("time_transform", "propose_ei_coal_tree_tt_rejectnegs.jl"))
include(srcdir("time_transform", "propose_ei_coal_tree_tt_rejectbadpop.jl"))

### First, simulate the EI MJP process ###
individ_results = DataFrame()
state_results = DataFrame()
I_init = 1
r0 = 2 # old one was 2
nu = 1/3
test_alpha = r0 * nu
gamma = 1/2
stop_time = 35.0
individ_frame = DataFrame()
state_frame = DataFrame()
candidate_seed = 235
total_pop_sims = 1000
all_individ_dfs = []
all_state_dfs = []
i = 1
while i <= total_pop_sims
    individ_frame, state_frame = sim_ei(I_init, test_alpha, nu, gamma, stop_time)
    # check to see if we have enough individuals to sample
    # in the isochronous case, max is 100, to be safe, assume we need 130
    if state_frame.I[end] >= 5
        individ_frame.sim .= i
        individ_frame.real_seed .= candidate_seed
        state_frame.sim .= i
        state_frame.real_seed .= candidate_seed
        push!(all_individ_dfs, individ_frame)
        push!(all_state_dfs, state_frame)
        i+=1
    end
    candidate_seed += 1
end
final_individ_frame = vcat(all_individ_dfs...)
final_state_frame = vcat(all_state_dfs...)

CSV.write(datadir("compare_data", "compare_ei_state_frame.csv"), final_state_frame)
CSV.write(datadir("compare_data", "compare_ei_individ_frame.csv"), final_individ_frame)

### now simulate 100 coalescent trees using deterministic coalescent construction ###
last_samp_time = stop_time
num_samps = 5
results5 = Array{Float64}(undef, 1000, 5) 
forward_samp_times = []
forward_samp_lin = []
forward_samp_ids = []
list_of_ids = []
for s in 1:1000
    current_individ_frame = final_individ_frame[final_individ_frame.sim .== s, :]
    individ_frame_clean = dropmissing(current_individ_frame, [:infectious_time])
    eligible_individuals = filter(row -> row.infectious_time < last_samp_time && (ismissing(row.recov_time) || row.recov_time > last_samp_time), eachrow(individ_frame_clean))
    samp_ids = sample(eligible_individuals.id, num_samps, replace = false)
    tree, tree_coal_pairs = create_ei_coal_tree_hetchron2(samp_ids, last_samp_time, forward_samp_times, 
    forward_samp_ids, current_individ_frame)
    coal_times = abs.(tree.time[tree.event .== "coal"] .- last_samp_time)
    # Create a new row with the iteration number and the results
    new_row = [s; coal_times...]
    # Add the new row to the results array
    results5[s, :] = new_row
end

# save the results
column_names = [:iteration, :coal_time1, :coal_time2, :coal_time3, :coal_time4]

# Convert the array to a DataFrame
results_df5 = DataFrame(results5, column_names)
results_df5.diff1 = results_df5.coal_time1
results_df5.diff2 = results_df5.coal_time2 .- results_df5.coal_time1
results_df5.diff3 = results_df5.coal_time3 .- results_df5.coal_time2
results_df5.diff4 = results_df5.coal_time4 .- results_df5.coal_time3

# Save the results
CSV.write(datadir("compare_data", "compare_ei_impcoal_results.csv"), results_df5)

### time transformation, rejecting negative coalescent times, resetting to last valid coal time ###
final_state_frame = CSV.read(datadir("compare_data", "compare_ei_state_frame.csv"), DataFrame)
last_samp_time = stop_time
num_samps = 5
# Add a new column reverse_time with the absolute difference from samp_time
init_nI = num_samps
init_nE = 0
max_time = last_samp_time
results5_noneg = Array{Float64}(undef, 1000, 5)  # 10000 rows, 5 columns

for s in 1:1000
    Random.seed!(s)
    print("s", s)
    comp_frame = filter(row -> row.time <= last_samp_time && row.sim == s, final_state_frame)
    comp_frame.reverse_time = abs.(comp_frame.time .- last_samp_time)
    reverse_times = sort(comp_frame.reverse_time)
    coal_times = propose_ei_coal_tree_tt_rejectbadpop(last_samp_time; samp_time = last_samp_time, reverse_times = reverse_times, 
    comp_traj = comp_frame, gamma = gamma, 
    alpha = test_alpha, init_nE = init_nE, init_nI = init_nI)[5]
    # Create a new row with the iteration number and the results
    new_row = [s; coal_times...]
    # Add the new row to the results array
    results5_noneg[s, :] = new_row
end
column_names = [:iteration, :coal_time1, :coal_time2, :coal_time3, :coal_time4]

results5_df5_ttnoneg = DataFrame(results5_noneg, column_names)
results5_df5_ttnoneg.diff1 = results5_df5_ttnoneg.coal_time1
results5_df5_ttnoneg.diff2 = results5_df5_ttnoneg.coal_time2 .- results5_df5_ttnoneg.coal_time1
results5_df5_ttnoneg.diff3 = results5_df5_ttnoneg.coal_time3 .- results5_df5_ttnoneg.coal_time2
results5_df5_ttnoneg.diff4 = results5_df5_ttnoneg.coal_time4 .- results5_df5_ttnoneg.coal_time3
CSV.write(datadir("compare_data", "compare_ei_tt_badpop_results.csv"), results5_df5_ttnoneg)

### time transformation, rejecting negative coalescent times, use ODE result instead of true traj ###
include(srcdir("new_thin", "sim_next_eicoal_time_tt.jl"))
include(srcdir("new_thin", "propose_ei_coal_tree_tt_rejectnegs.jl"))
final_state_frame = CSV.read(datadir("compare_data", "compare_ei_state_frame.csv"), DataFrame)
last_samp_time = stop_time
num_samps = 5
# Add a new column reverse_time with the absolute difference from samp_time
init_nI = num_samps
init_nE = 0
max_time = last_samp_time
results5_ode = Array{Float64}(undef, 1000, 5)  # 10000 rows, 5 columns
# calculates the solutions to the EI ode system
include(srcdir("calc_ei_trajectoriesv2.jl"))
# needed to calcluate the solutions to the EI ode system
include(srcdir("core_ei_func.jl"))
comp_times = vcat(collect(0.5:0.5:35))

E_traj = zeros(length(comp_times))
I_traj = zeros(length(comp_times))
init_e0 = 0
init_i0 = 1
calc_ei_trajectoriesv2!(comp_times, [0], [test_alpha], gamma, nu, init_e0,init_i0, E_traj, I_traj)
time = vcat(0, comp_times)
ode_comp_frame = DataFrame((time = time, E = vcat(init_e0, E_traj),I = vcat(init_i0, I_traj)))
ode_comp_frame.reverse_time = abs.(ode_comp_frame.time .- last_samp_time)
reverse_times = sort(ode_comp_frame.reverse_time)
for s in 1:1000
    Random.seed!(s)
    print("s", s)
    coal_times = propose_ei_coal_tree_tt_rejectnegs(last_samp_time; samp_time = last_samp_time, reverse_times = reverse_times, 
    comp_traj = ode_comp_frame, gamma = gamma, 
    alpha = test_alpha, init_nE = init_nE, init_nI = init_nI)[5]
    # Create a new row with the iteration number and the results
    new_row = [s; coal_times...]
    # Add the new row to the results array
    results5_ode[s, :] = new_row
end
column_names = [:iteration, :coal_time1, :coal_time2, :coal_time3, :coal_time4]

results_df5_ode = DataFrame(results5_ode, column_names)
results_df5_ode.diff1 = results_df5_ode.coal_time1
results_df5_ode.diff2 = results_df5_ode.coal_time2 .- results_df5_ode.coal_time1
results_df5_ode.diff3 = results_df5_ode.coal_time3 .- results_df5_ode.coal_time2
results_df5_ode.diff4 = results_df5_ode.coal_time4 .- results_df5_ode.coal_time3
CSV.write(datadir("compare_data", "compare_ei_tt_ode_results.csv"), results_df5_ode)


### time transformation, ODE, reject neg, use grids of each individ sim ###
include(srcdir("time_transform", "sim_next_eicoal_time_tt.jl"))
include(srcdir("time_transform", "propose_ei_coal_tree_tt_rejectnegs.jl"))
final_state_frame = CSV.read(datadir("compare_data", "compare_ei_state_frame.csv"), DataFrame)
last_samp_time = stop_time
num_samps = 5
# Add a new column reverse_time with the absolute difference from samp_time
init_nI = num_samps
init_nE = 0
max_time = last_samp_time
results5_ode = Array{Float64}(undef, 1000, 5)  # 10000 rows, 5 columns
# calculates the solutions to the EI ode system
include(srcdir("calc_ei_trajectoriesv2.jl"))
# needed to calcluate the solutions to the EI ode system
include(srcdir("core_ei_func.jl"))
comp_times = vcat(collect(0.5:0.5:35))

E_traj = zeros(length(comp_times))
I_traj = zeros(length(comp_times))
init_e0 = 0
init_i0 = 1
calc_ei_trajectoriesv2!(comp_times, [0], [test_alpha], gamma, nu, init_e0,init_i0, E_traj, I_traj)
time = vcat(0, comp_times)
ode_comp_frame = DataFrame((time = time, E = vcat(init_e0, E_traj),I = vcat(init_i0, I_traj)))
ode_comp_frame.reverse_time = abs.(ode_comp_frame.time .- last_samp_time)
reverse_times = sort(ode_comp_frame.reverse_time)
for s in 1:1000
    Random.seed!(s)
    print("s", s)
    comp_frame = filter(row -> row.time <= last_samp_time && row.sim == s, final_state_frame)
    comp_frame.reverse_time = abs.(comp_frame.time .- last_samp_time)
    reverse_times = sort(comp_frame.reverse_time)
    E_traj = zeros(length(comp_times))
    I_traj = zeros(length(comp_times))
    init_e0 = 0
    init_i0 = 1
    calc_ei_trajectoriesv2!(comp_times, [0], [test_alpha], gamma, nu, init_e0,init_i0, E_traj, I_traj)
    time = vcat(0, comp_times)
    ode_comp_frame = DataFrame((time = time, E = vcat(init_e0, E_traj),I = vcat(init_i0, I_traj)))

    coal_times = propose_ei_coal_tree_tt_rejectnegs(last_samp_time; samp_time = last_samp_time, reverse_times = reverse_times, 
    comp_traj = ode_comp_frame, gamma = gamma, 
    alpha = test_alpha, init_nE = init_nE, init_nI = init_nI)[5]
    # Create a new row with the iteration number and the results
    new_row = [s; coal_times...]
    # Add the new row to the results array
    results5_ode[s, :] = new_row
end
column_names = [:iteration, :coal_time1, :coal_time2, :coal_time3, :coal_time4]

results_df5_ode = DataFrame(results5_ode, column_names)
results_df5_ode.diff1 = results_df5_ode.coal_time1
results_df5_ode.diff2 = results_df5_ode.coal_time2 .- results_df5_ode.coal_time1
results_df5_ode.diff3 = results_df5_ode.coal_time3 .- results_df5_ode.coal_time2
results_df5_ode.diff4 = results_df5_ode.coal_time4 .- results_df5_ode.coal_time3
CSV.write(datadir("compare_data", "compare_ei_tt_ode_results.csv"), results_df5_ode)


