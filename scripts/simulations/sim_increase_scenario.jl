using DrWatson
using CSV
using Distributions
using Statistics
using StatsBase
using Combinatorics
using Random
include(srcdir("sim_ei.jl"))
include(srcdir("sim_seir.jl"))
include(srcdir("create_coal_trees.jl"))
include(srcdir("construct_newick_tree.jl"))

### how many sims? ###
num_pop_sims = 1
### Simulate an SEIR trajectory ###
N = 15000
I_init = 1
r0_vec = vcat(exp.(collect(log(1.3):(log(2) - log(1.3))/(5*7):log(2))), 2)
nu = 1/7
beta_vec = r0_vec .* nu
beta_times = collect(0.0:1.0:length(r0_vec)-1)
gamma = 1/4
stop_time = 22.0 * 7
last_samp_time = stop_time - 1
i = 0
seed_num = 100
all_individ_dfs = []
all_state_dfs = []
while i < num_pop_sims
    Random.seed!(seed_num)
    individ_frame, state_frame = sim_timevarying_seir(N, I_init, beta_vec, beta_times, nu, gamma, stop_time)

    # check to see if we have enough individuals to sample
    # in the isochronous case, max is 100, to be safe, assume we need 120
    # find the index closes to last_samp_time
    if maximum(state_frame.time) >= stop_time
        last_samp_index = findfirst(x -> x >= last_samp_time, state_frame.time)
    else
        last_samp_index = size(state_frame, 1)
    end
    if state_frame.I[last_samp_index] >= 130
        i += 1
        individ_frame.sim .= i
        individ_frame.real_seed .= seed_num+1
        state_frame.sim .= i
        state_frame.real_seed .= seed_num+1
        push!(all_individ_dfs, individ_frame)
        push!(all_state_dfs, state_frame)
    end
    seed_num += 1
end
# save the individ frame and the state frame
final_individ_frame = vcat(all_individ_dfs...)
final_state_frame = vcat(all_state_dfs...)
CSV.write(datadir("sim_data", "increase_individ_frame.csv"), final_individ_frame)
CSV.write(datadir("sim_data", "increase_state_frame.csv"), final_state_frame)
### sampling scheme 1 isochronous 50 ###
num_samp_sims = 50
num_samps = 50
iso = true 
forward_samp_times = []
forward_samp_lin = []
forward_samp_ids = []
all_iso50_trees = []
individ_frame = CSV.read(datadir("sim_data", "increase_individ_frame.csv"), DataFrame)
individ_frame.history = [eval(Meta.parse(individ_frame.history[i])) for i in 1:size(individ_frame, 1)]

for j in 1:num_samp_sims
    # sample infectious individuals from time last_samp_time 
    Random.seed!(j)
    last_samp_ids = sample_infected_ids(individ_frame, last_samp_time, forward_samp_times, forward_samp_lin, num_samps, iso)[1]
    # create a tree from the samp_ids
    tree, tree_coal_pairs = create_ei_coal_tree_hetchron2(last_samp_ids, last_samp_time, forward_samp_times, forward_samp_ids, individ_frame)
    # record the sim and seed number
    tree.sim .= j
    tree.real_seed .= j
    tree.time = tree.time = abs.(tree.time .- last_samp_time)

    # store the tree 
    push!(all_iso50_trees, tree)
end
final_iso50_tree = vcat(all_iso50_trees...)
# save the trees 
CSV.write(datadir("sim_data", "increase_iso50_trees.csv"), final_iso50_tree)

### sampling scheme 2 isochronous 100 ###
num_samp_sims = 50
num_samps = 100
iso = true 
forward_samp_times = []
forward_samp_lin = []
forward_samp_ids = []
all_iso100_trees = []
individ_frame = CSV.read(datadir("sim_data", "increase_individ_frame.csv"), DataFrame)
individ_frame.history = [eval(Meta.parse(individ_frame.history[i])) for i in 1:size(individ_frame, 1)]

for j in 1:num_samp_sims
    # sample infectious individuals from time last_samp_time 
    Random.seed!(j)
    last_samp_ids = sample_infected_ids(individ_frame, last_samp_time, forward_samp_times, forward_samp_lin, num_samps, iso)[1]
    # create a tree from the samp_ids
    tree, tree_coal_pairs = create_ei_coal_tree_hetchron2(last_samp_ids, last_samp_time, forward_samp_times, forward_samp_ids, individ_frame)
    # record the sim and seed number
    tree.sim .= j
    tree.real_seed .= j
    tree.time = tree.time = abs.(tree.time .- last_samp_time)

    # store the tree 
    push!(all_iso100_trees, tree)
end
final_iso100_tree = vcat(all_iso100_trees...)
# save the trees 
CSV.write(datadir("sim_data", "increase_iso100_trees.csv"), final_iso100_tree)

### sampling scheme 3 heterochronous 50 ###
num_samp_sims = 50
total_samps = 50
last_num_samps = 30
individ_frame = CSV.read(datadir("sim_data", "increase_individ_frame.csv"), DataFrame)
individ_frame.history = [eval(Meta.parse(individ_frame.history[i])) for i in 1:size(individ_frame, 1)]

# five weeks before last sampling time, we allow sampling
samp_window_start = last_samp_time - 5*7
iso = false 
forward_samp_times = []
forward_samp_lin = ones(Int, total_samps - last_num_samps)
forward_samp_ids = []
all_het50_trees = []
# individ_frame = all_individ_dfs[1]
for j in 1:num_samp_sims
    # sample infectious individuals from time last_samp_time 
    Random.seed!(j)
    # sample the times at which sampling occurs
    forward_samp_times = sort(rand(Uniform(samp_window_start, last_samp_time), total_samps - last_num_samps))
    # sample the ids of those who were sampled at times forward_samp_times and last_samp_time
    samp_ids = sample_infected_ids(individ_frame, last_samp_time, forward_samp_times, forward_samp_lin, last_num_samps, iso)
    last_samp_ids = samp_ids[end]
    forward_samp_ids = samp_ids[1:end-1]
    # create a tree from the samp_ids
    tree, tree_coal_pairs = create_ei_coal_tree_hetchron2(last_samp_ids, last_samp_time, forward_samp_times, forward_samp_ids, individ_frame)
    # record the sim and seed number
    tree.sim .= j
    tree.real_seed .= j
    tree.time = tree.time = abs.(tree.time .- last_samp_time)

    # store the tree 
    push!(all_het50_trees, tree)
end
final_het50_tree = vcat(all_het50_trees...)
# save the trees
CSV.write(datadir("sim_data", "increase_het50_trees.csv"), final_het50_tree)

### sampling scheme 4 heterochronous 100 ###
num_samp_sims = 50
total_samps = 100
last_num_samps = 70
# five weeks before last sampling time, we allow sampling
samp_window_start = last_samp_time - 5*7
iso = false 
forward_samp_times = []
forward_samp_lin = ones(Int, total_samps - last_num_samps)
forward_samp_ids = []
all_het100_trees = []
individ_frame = CSV.read(datadir("sim_data", "increase_individ_frame.csv"), DataFrame)
individ_frame.history = [eval(Meta.parse(individ_frame.history[i])) for i in 1:size(individ_frame, 1)]
for j in 1:num_samp_sims
    # sample infectious individuals from time last_samp_time 
    Random.seed!(j)
    # sample the times at which sampling occurs
    forward_samp_times = sort(rand(Uniform(samp_window_start, last_samp_time), total_samps - last_num_samps))
    # sample the ids of those who were sampled at times forward_samp_times and last_samp_time
    samp_ids = sample_infected_ids(individ_frame, last_samp_time, forward_samp_times, forward_samp_lin, last_num_samps, iso)
    last_samp_ids = samp_ids[end]
    forward_samp_ids = samp_ids[1:end-1]
    # create a tree from the samp_ids
    tree, tree_coal_pairs = create_ei_coal_tree_hetchron2(last_samp_ids, last_samp_time, forward_samp_times, forward_samp_ids, individ_frame)
    # record the sim and seed number
    tree.sim .= j
    tree.real_seed .= j
    tree.time = tree.time = abs.(tree.time .- last_samp_time)

    # store the tree 
    push!(all_het100_trees, tree)
end
final_het100_tree = vcat(all_het100_trees...)
# save the trees
CSV.write(datadir("sim_data", "increase_het100_trees.csv"), final_het100_tree)

### update sim_dict.csv ###
sim_dict = CSV.read(datadir("sim_data", "sim_dict.csv"), DataFrame)
# filter old versions of sim_id 5,6,7,8
sim_dict = sim_dict[.!(sim_dict.sim_id .∈ Ref([5, 6, 7, 8])), :]
new_sim_dict = DataFrame(
    sim_id = [5,6,7,8],
    sim = ["increase_iso50", "increase_iso100", "increase_het50", "increase_het100"],
    individ_frame_file = ["increase_individ_frame.csv", "increase_individ_frame.csv", "increase_individ_frame.csv", "increase_individ_frame.csv"],
    state_frame_file = ["increase_state_frame.csv", "increase_state_frame.csv", "increase_state_frame.csv", "increase_state_frame.csv"],
    tree_file = ["increase_iso50_trees.csv", "increase_iso100_trees.csv", "increase_het50_trees.csv", "increase_het100_trees.csv"],
    sim_julia_file = ["sim_increase_scenario.jl", "sim_increase_scenario.jl", "sim_increase_scenario.jl", "sim_increase_scenario.jl"],
    N = [N, N, N, N],
    I_init = [1, 1, 1, 1],
    r0 = [r0_vec[1], r0_vec[1], r0_vec[1], r0_vec[1]],
    nu = [nu, nu, nu, nu],
    gamma = [gamma, gamma, gamma, gamma],
    stop_time = [22.0 * 7, 22.0 * 7, 22.0 * 7, 22.0 * 7],
    last_samp_time = [22.0 * 7 - 1, 22.0 * 7 - 1, 22.0 * 7 - 1, 22.0 * 7 - 1]
)
sim_dict = vcat(sim_dict, new_sim_dict)
CSV.write(datadir("sim_data", "sim_dict.csv"), sim_dict)

### make newick tree files
# sim 5
trees = CSV.read(datadir("sim_data", "increase_iso50_trees.csv"), DataFrame)
trees.og_ids = [eval(Meta.parse(trees.og_ids[i])) for i in 1:size(trees, 1)]
# filter for sim == 1
coal_and_samp_tree = trees[trees.sim .== 1, :]
reverse_samp_lin = [50] 
newick_tree = construct_newick_tree(coal_and_samp_tree, reverse_samp_lin)
# save the newick tree as .tree
write(datadir("sim_data", "increase_iso50_simnum1.tree"), newick_tree)
# sim 6
trees = CSV.read(datadir("sim_data", "increase_iso100_trees.csv"), DataFrame)
trees.og_ids = [eval(Meta.parse(trees.og_ids[i])) for i in 1:size(trees, 1)]
# filter for sim == 1
coal_and_samp_tree = trees[trees.sim .== 1, :]
reverse_samp_lin = [100]
newick_tree = construct_newick_tree(coal_and_samp_tree, reverse_samp_lin)
# save the newick tree as .tree
write(datadir("sim_data", "increase_iso100_simnum1.tree"), newick_tree)
# sim 7
trees = CSV.read(datadir("sim_data", "increase_het50_trees.csv"), DataFrame)
trees.og_ids = [eval(Meta.parse(trees.og_ids[i])) for i in 1:size(trees, 1)]
# filter for sim == 1
coal_and_samp_tree = trees[trees.sim .== 1, :]
reverse_samp_lin = vcat(30, ones(Int, 20))
newick_tree = construct_newick_tree(coal_and_samp_tree, reverse_samp_lin)
# save the newick tree as .tree
write(datadir("sim_data", "increase_het50_simnum1.tree"), newick_tree)
# sim 8
trees = CSV.read(datadir("sim_data", "increase_het100_trees.csv"), DataFrame)
trees.og_ids = [eval(Meta.parse(trees.og_ids[i])) for i in 1:size(trees, 1)]
# filter for sim == 1
coal_and_samp_tree = trees[trees.sim .== 1, :]
reverse_samp_lin = vcat(70, ones(Int, 30))
newick_tree = construct_newick_tree(coal_and_samp_tree, reverse_samp_lin)
# save the newick tree as .tree
write(datadir("sim_data", "increase_het100_simnum1.tree"), newick_tree)
# sim 7 simnum2
trees = CSV.read(datadir("sim_data", "increase_het50_trees.csv"), DataFrame)
trees.og_ids = [eval(Meta.parse(trees.og_ids[i])) for i in 1:size(trees, 1)]
# filter for sim == 1
coal_and_samp_tree = trees[trees.sim .== 2, :]
reverse_samp_lin = vcat(30, ones(Int, 20))
newick_tree = construct_newick_tree(coal_and_samp_tree, reverse_samp_lin)
# save the newick tree as .tree
write(datadir("sim_data", "increase_het50_simnum2.tree"), newick_tree)
# sim 5 simnum 10
trees = CSV.read(datadir("sim_data", "increase_iso50_trees.csv"), DataFrame)
trees.og_ids = [eval(Meta.parse(trees.og_ids[i])) for i in 1:size(trees, 1)]
# filter for sim == 1
coal_and_samp_tree = trees[trees.sim .== 10, :]
reverse_samp_lin = [50]
newick_tree = construct_newick_tree(coal_and_samp_tree, reverse_samp_lin)
# save the newick tree as .tree
write(datadir("sim_data", "increase_iso50_simnum10.tree"), newick_tree)
# sim 5 simnum 2
trees = CSV.read(datadir("sim_data", "increase_iso50_trees.csv"), DataFrame)
trees.og_ids = [eval(Meta.parse(trees.og_ids[i])) for i in 1:size(trees, 1)]
# filter for sim == 2
coal_and_samp_tree = trees[trees.sim .== 2, :]
reverse_samp_lin = [50]
newick_tree = construct_newick_tree(coal_and_samp_tree, reverse_samp_lin)
# save the newick tree as .tree
write(datadir("sim_data", "increase_iso50_simnum2.tree"), newick_tree)
# sim5 simnum 3
trees = CSV.read(datadir("sim_data", "increase_iso50_trees.csv"), DataFrame)
trees.og_ids = [eval(Meta.parse(trees.og_ids[i])) for i in 1:size(trees, 1)]
# filter for sim == 3
coal_and_samp_tree = trees[trees.sim .== 3, :]
reverse_samp_lin = [50]
newick_tree = construct_newick_tree(coal_and_samp_tree, reverse_samp_lin)
# save the newick tree as .tree
write(datadir("sim_data", "increase_iso50_simnum3.tree"), newick_tree)
