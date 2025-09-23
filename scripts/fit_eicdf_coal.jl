# fit partially observed EI coalescent model to simulated data
using CSV
using DataFrames
using DrWatson
using Distributions
using Revise
using JLD2
using Random
using LinearAlgebra
using ExponentialUtilities
using StatsBase
# calculates the solutions to the EI ode system
include(srcdir("calc_ei_trajectoriesv2.jl"))
# needed to calcluate the solutions to the EI ode system
include(srcdir("core_ei_func.jl"))
# node likelihood functions
include(srcdir("node_simple", "sample_init_params_simple.jl"))
include(srcdir("node_simple", "sample_nodescdf_andparams.jl"))
include(srcdir("node_simple", "sample_ess_simple.jl"))
include(srcdir("node_simple", "calc_node_loglik_simplev2.jl"))
include(srcdir("node_simple", "calc_node_loglik_simple_safev2.jl"))
include(srcdir("node_simple", "sample_internal_nodesv2.jl"))
num_samples = 350000
discard_initial = Int(round(num_samples/2))
num_thin = 35
# decide on the simulation id and sim number 
sim_id =
if length(ARGS) == 0
    6
else
  parse(Int64, ARGS[1])
end

sim_num = 
if length(ARGS) == 0
  4
else 
  parse(Int64, ARGS[2])
end 
# load the simulated data and priors 
include(projectdir("src/load_dp_ei_coal.jl"))
# choose alpha_times 
change_grid = 7.0
samp_time = maximum(tree.time) 
if samp_time % change_grid == 0
    param_change_max = maximum(samp_time) - change_grid
  else 
    param_change_max = maximum(samp_time)
  end 
alpha_times = vcat(0,collect(change_grid:change_grid:param_change_max))
# find the coal and sample times
coal_and_samp_tree = tree[(tree.event .== "coal") .| (tree.event .== "samp"),:]
coal_times = coal_and_samp_tree.time[coal_and_samp_tree.event .== "coal"]
# remove the first sampling time, that's just the start time
reverse_samp_times = coal_and_samp_tree.time[coal_and_samp_tree.event .== "samp"][2:end]
coal_and_samp_times = sort(vcat(reverse_samp_times, coal_times))
coal_times = sort!(coal_times)
reverse_samp_times = sort!(reverse_samp_times)

# choose comp_times to be the alpha times and the times of coal and sample events
comp_times = sort(union(alpha_times[2:end], abs.(coal_times[1:end-1] .- samp_time), abs.(reverse_samp_times .- samp_time), samp_time))
reverse_times = reverse(abs.(vcat(0.0, comp_times) .- samp_time))
# est_times is the times at which we will estimate the states 
est_times = sort(union(coal_times, reverse_samp_times, abs.(alpha_times[2:end] .- samp_time)))
start_time = 0.0
num_lineages = tree.nE[1] + tree.nI[1]
reverse_alpha_times = reverse(abs.(alpha_times[2:end] .- samp_time))
curr_lin = create_curr_lin_vec(num_lineages, coal_times, reverse_samp_times, reverse_samp_lin, reverse_alpha_times)

# choose initial values for the states and the parameters
Random.seed!(sim_num)
init_rw_sigma, init_rts, init_gamma, init_nu, init_e0, init_i0 = sample_initial_params_simple(log_rw_mean, log_rw_sigma_sd, log_rt_init_mean, log_rt_init_sd,
    log_gamma_mean, log_gamma_sd, log_nu_mean, log_nu_sd, log_e0_mean, log_e0_sd, log_i0_mean, log_i0_sd, 
  alpha_times, curr_lin)
alpha_vec = init_rts .* init_nu
E_traj = zeros(length(comp_times))
I_traj = zeros(length(comp_times))
calc_ei_trajectoriesv2!(comp_times, alpha_times, alpha_vec, init_gamma, init_nu, init_e0,init_i0, E_traj, I_traj)
reverse_E = vcat(reverse(E_traj), init_e0)
reverse_I = vcat(reverse(I_traj), init_i0)
total_pop = reverse_E .+ reverse_I
max_lineages = Int(maximum(curr_lin))
init_dist = zeros(max_lineages + 1)
mat_size = 50
last_samp_time = samp_time
# this is the cutoff for using the krylov method, if the time step is smaller than this, we use the krylov method
# for the even simulations, sometimes they fail with this large cutoff, we will use 0.5 instead
tstep_cutoff = 0.5 
Random.seed!(sim_num)
log_lik, log_lik_vec, new_sampled_node_states =sample_internal_nodes_safev2(num_lineages, est_times, coal_times, 
init_dist, start_time, last_samp_time, reverse_samp_times, reverse_samp_lin, alpha_times, init_gamma, alpha_vec, 
reverse_E, reverse_I, max_lineages, mat_size)

# fit the model
est_states = new_sampled_node_states
natural_vars = vcat(init_gamma, init_nu, init_e0, init_i0,  init_rw_sigma, init_rts[1], init_rts[2:end])
q_cur = log.(natural_vars) .- vcat(log_prior_means, repeat([log_rt_init_mean], length(init_rts)-1))
l_cur = log_lik
Random.seed!(sim_num)
@time my_samples, my_states = sample_nodescdf_andparams!(q_cur, l_cur, cholC, log_prior_means, num_lineages, est_times, coal_times, 
    est_states, start_time,last_samp_time, reverse_samp_times, reverse_samp_lin, alpha_times, mat_size, curr_lin,
     num_samples, discard_initial, num_thin, tstep_cutoff)
# make dataframe of results
rt_columns = ["rt_t_values[$i]" for i in 0:(length(alpha_times)-1)]
other_columns = [:gamma, :nu, :e0, :i0, :rw_sigma]
# rewrite other_columns as characters
other_columns = [string(col) for col in other_columns]
all_columns = vcat(other_columns, rt_columns, ["log_likelihood", "actual_iteration"])
my_samples_frame = DataFrame(my_samples, all_columns)
# save results
CSV.write(resultsdir("my_generated_quantities", 
string("ei_cdf_sim", sim_id, "_simnum", sim_num, ".csv")), my_samples_frame)
