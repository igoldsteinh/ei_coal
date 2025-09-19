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
using Optim
using ExponentialUtilities
include(srcdir("sim_ei.jl"))
include(srcdir("create_coal_trees.jl"))
include(srcdir("new_thin", "propose_ei_coal_tree.jl"))
include(srcdir("new_thin", "sim_ei_coal_tree_thin.jl"))
include(srcdir("new_thin", "ei_event_rate_sampler.jl"))
include(srcdir("new_thin", "calc_ei_probs_julia.jl"))
include(srcdir("new_thin", "sim_next_eicoal_time_tt.jl"))
include(srcdir("new_thin", "propose_ei_coal_tree_tt.jl"))
include(srcdir("new_thin", "sim_next_eicoal_time_tt_mod.jl"))
include(srcdir("new_thin", "propose_ei_coal_tree_tt_mod.jl"))
include(srcdir("new_thin", "propose_ei_coal_tree_tt_naive.jl"))
### First, simulate the EI MJP process ###
individ_results = DataFrame()
state_results = DataFrame()

### Baseline simulation ###
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
# Random.seed!(candidate_seed)
# print(candidate_seed)
# individ_frame, state_frame = sim_ei(I_init, test_alpha, nu, gamma, stop_time)
# maximum(state_frame.E .+ state_frame.I) 
# maximum(individ_frame.id)
final_individ_frame = vcat(all_individ_dfs...)
final_state_frame = vcat(all_state_dfs...)

CSV.write(datadir("compare_data", "compare_ei_state_frame.csv"), final_state_frame)
CSV.write(datadir("compare_data", "compare_ei_individ_frame.csv"), final_individ_frame)

### now simulate 100 coalescent trees using deterministic coalescent construction ###
# final_individ_frame = CSV.read(datadir("compare_data", "compare_ei_individ_frame.csv"), DataFrame)
# final_individ_frame.history = [eval(Meta.parse(final_individ_frame.history[i])) for i in 1:size(final_individ_frame, 1)]

last_samp_time = stop_time
num_samps = 5
results5 = Array{Float64}(undef, 1000, 5)  # 1000 rows, 5 columns# simulate many coalescences for 5 individuals
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

### Save the results
CSV.write(datadir("compare_data", "compare_ei_impcoal_results.csv"), results_df5)


### do the same but with time transformation ###
final_state_frame = CSV.read(datadir("compare_data", "compare_ei_state_frame.csv"), DataFrame)
last_samp_time = stop_time
num_samps = 5
# Add a new column reverse_time with the absolute difference from samp_time
init_nI = num_samps
init_nE = 0
max_time = last_samp_time
results5_tt = Array{Float64}(undef, 1000, 5)  # 10000 rows, 5 columns

for s in 1:1000
    Random.seed!(s)
    print("s", s)
    comp_frame = filter(row -> row.time <= last_samp_time && row.sim == s, final_state_frame)
    comp_frame.reverse_time = abs.(comp_frame.time .- last_samp_time)
    reverse_times = sort(comp_frame.reverse_time)
    coal_times = propose_ei_coal_tree_tt(last_samp_time; samp_time = last_samp_time, reverse_times = reverse_times, 
    comp_traj = comp_frame, gamma = gamma, 
    alpha = test_alpha, init_nE = init_nE, init_nI = init_nI)[7]
    # Create a new row with the iteration number and the results
    new_row = [s; coal_times...]
    # Add the new row to the results array
    results5_tt[s, :] = new_row
end
column_names = [:iteration, :coal_time1, :coal_time2, :coal_time3, :coal_time4]

results5_df5_tt = DataFrame(results5_tt, column_names)
results5_df5_tt.diff1 = results5_df5_tt.coal_time1
results5_df5_tt.diff2 = results5_df5_tt.coal_time2 .- results5_df5_tt.coal_time1
results5_df5_tt.diff3 = results5_df5_tt.coal_time3 .- results5_df5_tt.coal_time2
results5_df5_tt.diff4 = results5_df5_tt.coal_time4 .- results5_df5_tt.coal_time3
CSV.write(datadir("compare_data", "compare_ei_tt_results.csv"), results5_df5_tt)

# mean(results5_df5_tt.diff1)
# mean(results5_df5_tt.diff2)
# mean(results5_df5_tt.diff3)
# mean(results5_df5_tt.diff4)

### now redo time transformation with naive rejection sampling ###
include(srcdir("new_thin", "sim_next_eicoal_time_tt_naive.jl"))
include(srcdir("new_thin", "propose_ei_coal_tree_tt_naive3.jl"))
final_state_frame = CSV.read(datadir("compare_data", "compare_ei_state_frame.csv"), DataFrame)
last_samp_time = stop_time
num_samps = 5
# Add a new column reverse_time with the absolute difference from samp_time
init_nI = num_samps
init_nE = 0
max_time = last_samp_time
results5_ttnaive = Array{Float64}(undef, 1000, 5)  # 10000 rows, 5 columns

for s in 1:1000
    Random.seed!(s)
    print("s", s)
    comp_frame = filter(row -> row.time <= last_samp_time && row.sim == s, final_state_frame)
    comp_frame.reverse_time = abs.(comp_frame.time .- last_samp_time)
    reverse_times = sort(comp_frame.reverse_time)
    coal_times = propose_ei_coal_tree_tt_naive3(last_samp_time; samp_time = last_samp_time, reverse_times = reverse_times, 
    comp_traj = comp_frame, gamma = gamma, 
    alpha = test_alpha, init_nE = init_nE, init_nI = init_nI)[5]
    # Create a new row with the iteration number and the results
    new_row = [s; coal_times...]
    # Add the new row to the results array
    results5_ttnaive[s, :] = new_row
end
column_names = [:iteration, :coal_time1, :coal_time2, :coal_time3, :coal_time4]

results5_df5_ttnaive = DataFrame(results5_ttnaive, column_names)
results5_df5_ttnaive.diff1 = results5_df5_ttnaive.coal_time1
results5_df5_ttnaive.diff2 = results5_df5_ttnaive.coal_time2 .- results5_df5_ttnaive.coal_time1
results5_df5_ttnaive.diff3 = results5_df5_ttnaive.coal_time3 .- results5_df5_ttnaive.coal_time2
results5_df5_ttnaive.diff4 = results5_df5_ttnaive.coal_time4 .- results5_df5_ttnaive.coal_time3
CSV.write(datadir("compare_data", "compare_ei_tt_naive_results.csv"), results5_df5_ttnaive)


### time transformation, rejecting negative coalescent times, resetting to last valid coal time ###
include(srcdir("new_thin", "sim_next_eicoal_time_tt.jl"))
include(srcdir("new_thin", "propose_ei_coal_tree_tt_rejectnegs.jl"))
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
    coal_times = propose_ei_coal_tree_tt_rejectnegs(last_samp_time; samp_time = last_samp_time, reverse_times = reverse_times, 
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
CSV.write(datadir("compare_data", "compare_ei_tt_noneg_results.csv"), results5_df5_ttnoneg)


# mean(results5_df5_tt_naive.diff1)
# mean(results5_df5_tt_naive.diff2)
# mean(results5_df5_tt_naive.diff3)
# mean(results5_df5_tt_naive.diff4)

### playing with Gillespie again ###
include(srcdir("new_thin", "sim_ei_coal_tree_gillespie.jl"))
state_frame = CSV.read(datadir("compare_data", "compare_ei_state_frame.csv"), DataFrame)
last_samp_time = stop_time
num_samps = 5

comp_frame = filter(row -> row.time <= last_samp_time, state_frame)
# Add a new column reverse_time with the absolute difference from samp_time
comp_frame.reverse_time = abs.(comp_frame.time .- last_samp_time)
reverse_times = sort(comp_frame.reverse_time)
init_nI = num_samps
init_nE = 0
max_time = last_samp_time
reverse_E = reverse(state_frame.E)
reverse_I = reverse(state_frame.I)
results5_gill = Array{Float64}(undef, 10000, 5)  # 10000 rows, 2 columns
for s in 1:10000
    Random.seed!(s)
    coal_times = sim_ei_coal_tree_gillespie(max_time; reverse_times, reverse_E, reverse_I, alpha, gamma, init_nE, init_nI)[5]
    new_row = [s; coal_times...]
    results5_gill[s, :] = new_row
end
column_names = [:iteration, :coal_time1]

results5_df5_gill = DataFrame(results5_gill, column_names)
results5_df5_gill.diff1 = results5_df5_gill.coal_time1
results5_df5_gill.diff2 = results5_df5_gill.coal_time2 .- results5_df5_gill.coal_time1
results5_df5_gill.diff3 = results5_df5_gill.coal_time3 .- results5_df5_gill.coal_time2
results5_df5_gill.diff4 = results5_df5_gill.coal_time4 .- results5_df5_gill.coal_time3
CSV.write(datadir("compare_data", "compare_ei_gill_results.csv"), results5_df5_gill)

# Current time=27.06691257935565
# Current nE= -1
# Current nI=4
# Current E=1.0
# Current I:=3.0
# next_time=27.357253377358127
# comp_idx=6128
# times=[0.0, 0.008711970732737342, 2.162918851614198, 2.58866813394578, 2.5958596879895874, 3.212059984687553, 3.4331142933696794, 4.781426673215429, 4.935873880121345, 5.05280952563978, 5.262618017462043, 5.7066384119906575, 6.3420674872004925, 6.402692821737721, 7.1555334401628325, 7.542075952616206, 7.982063941610932, 8.457959063901082, 9.40187503424005, 10.691631381739972, 11.050517646792388, 12.148427154689294, 12.596225504959744, 13.223056866841713, 14.332987693608809, 14.855609530199413, 14.938171997608013, 15.979840002743073, 16.024925575257104, 16.49866342066681, 16.542956541763296, 16.744691447229265, 16.751430218045495, 16.97848393411896, 17.5664608107382, 17.70338219543645, 17.799828869637413, 18.27808676643733, 19.418728513440918, 19.466280935301544, 20.130812063839013, 20.451381991181705, 21.026549406509695, 22.01036408731002, 22.019517722530846, 23.65170253378393, 23.668620129746877, 24.323795560656812, 25.017206743224484, 25.92435354200007, 26.053415463144333, 26.07019877483848]
# types=["start", "I2E", "I2E", "E2I", "I2E", "E2I", "E2I", "I2E", "I2E", "E2I", "I2E", "E2I", "I2E", "E2I", "E2I", "I2E", "E2I", "I2E", "E2I", "I2E", "E2I", "I2E", "I2E", "E2I", "E2I", "I2E", "I2E", "E2I", "I2E", "I2E", "E2I", "E2I", "I2E", "I2E", "I2E", "E2I", "E2I", "E2I", "I2E", "coal", "I2E", "I2E", "E2I", "I2E", "coal", "E2I", "E2I", "I2E", "I2E", "E2I", "E2I", "E2I"]
# nE_traj=[0, 1, 2, 1, 2, 1, 0, 1, 2, 1, 2, 1, 2, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 1, 0, 1, 2, 1, 2, 3, 2, 1, 2, 3, 4, 3, 2, 1, 2, 1, 2, 3, 2, 3, 2, 1, 0, 1, 2, 1, 0, -1]
# nI_traj=[5, 4, 3, 4, 3, 4, 5, 4, 3, 4, 3, 4, 3, 4, 5, 4, 5, 4, 5, 4, 5, 4, 3, 4, 5, 4, 3, 4, 3, 2, 3, 4, 3, 2, 1, 2, 3, 4, 3, 3, 2, 1, 2, 1, 1, 2, 3, 2, 1, 2, 3, 4]