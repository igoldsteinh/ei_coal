# file for converting simulated ei coal trees to newick format 
# for all 50 iso50 control trees used by BDMM
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
trees = CSV.read(datadir("sim_data", "control_iso50_trees.csv"), DataFrame)
trees.og_ids = [eval(Meta.parse(trees.og_ids[i])) for i in 1:size(trees, 1)]
# filter for sim == 1
for i in 2:50

    sim_val = i
    full_tree = trees[trees.sim .== sim_val, :]
    coal_and_samp_tree = full_tree[full_tree.event .== "samp" .|| full_tree.event .== "coal", :]
    reverse_samp_lin = vcat(50)
    newick_tree = construct_newick_tree(coal_and_samp_tree, reverse_samp_lin)
    # save the newick tree as .tree
    write(datadir("sim_data", string("control_iso50_simnum", i, ".tree")), newick_tree)
end
