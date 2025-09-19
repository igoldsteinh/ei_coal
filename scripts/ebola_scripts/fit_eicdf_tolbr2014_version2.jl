using CSV
using DataFrames
using DrWatson
using Distributions
using Revise
using JLD2
using Random
using LinearAlgebra
using MCMCChains
using ExponentialUtilities
using DynamicPPL
using StatsBase
using Plots

# helper functions for processing turing model
include(srcdir("phylo_ww.jl"))
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
# read in the data
samp_times = CSV.read(datadir("real_data", "ebola_lbrmcc_2014_samp_times.csv"), DataFrame)
coal_times_df = CSV.read(datadir("real_data", "ebola_lbrmcc_2014_coal_times.csv"), DataFrame)

coal_times = round.(coal_times_df.coal_times_days, digits = 4)
reverse_samp_times = samp_times.samp_times_4dig[2:end]
reverse_samp_lin = samp_times.samp_lin[2:end]
start_time = samp_times.samp_times_4dig[1]
num_lineages = samp_times.samp_lin[1]
last_samp_time = maximum(coal_times)
coal_and_samp_times = sort(vcat(reverse_samp_times, coal_times))
# pick weekly alpha change times
change_grid = 7.0
if last_samp_time % change_grid == 0
    param_change_max = last_samp_time - change_grid
  else 
    param_change_max = last_samp_time
  end 
alpha_times = vcat(0,collect(change_grid:change_grid:param_change_max))
coal_and_samp_times = sort(vcat(reverse_samp_times, coal_times))
coal_times = sort!(coal_times)
reverse_samp_times = sort!(reverse_samp_times)
samp_time = maximum(coal_times)
comp_times = sort(union(alpha_times[2:end], abs.(coal_times[1:end-1] .- samp_time), abs.(reverse_samp_times .- samp_time), samp_time))
reverse_times = reverse(abs.(vcat(0.0, comp_times) .- samp_time))
# est_times is the times at which we will estimate the states 
est_times = sort(union(coal_times, reverse_samp_times, abs.(alpha_times[2:end] .- samp_time)))
start_time = 0.0
reverse_alpha_times = reverse(abs.(alpha_times[2:end] .- samp_time))
curr_lin = create_curr_lin_vec(num_lineages, coal_times, reverse_samp_times, reverse_samp_lin, reverse_alpha_times)

### pick intial params ###
## set up priors
log_rt_init_mean = 0.7
log_rt_init_sd = 0.5
log_gamma_mean = log(1/7)
log_gamma_sd = 0.3
log_nu_mean = log(1/7)
log_nu_sd = 0.3
# note we're hardcoding in e0 has to be above 1
log_e0_mean = log(1.1)
log_e0_sd = 0.05
log_i0_mean = log(1.1)
log_i0_sd = 0.05
log_rw_mean = log(0.05)
log_rw_sigma_sd = 0.2

log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
Random.seed!(123467)
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
tstep_cutoff = 0.5 # this is the cutoff for using the krylov method, if the time step is larger than this, we don't use the krylov method
Random.seed!(123467)
log_lik, log_lik_vec, new_sampled_node_states =sample_internal_nodes_safev2(num_lineages, est_times, coal_times, 
init_dist, start_time, last_samp_time, reverse_samp_times, reverse_samp_lin, alpha_times, init_gamma, alpha_vec, 
reverse_E, reverse_I, max_lineages, mat_size)

est_states = new_sampled_node_states
natural_vars = vcat(init_gamma, init_nu, init_e0, init_i0,  init_rw_sigma, init_rts[1], init_rts[2:end])
q_cur = log.(natural_vars) .- vcat(log_prior_means, repeat([log_rt_init_mean], length(init_rts)-1))
l_cur = log_lik
num_samples = 900000 
num_thin = 90
discard_initial = 0
discard_initial = Int(round(num_samples/2))
Random.seed!(123467)
@time my_samples, my_states = sample_nodescdf_andparams!(q_cur, l_cur, cholC, log_prior_means, num_lineages, est_times, coal_times, 
    est_states, start_time,last_samp_time, reverse_samp_times, reverse_samp_lin, alpha_times, mat_size, curr_lin,
     num_samples, discard_initial, num_thin, tstep_cutoff)

# well lets see what it looks like I guess 
rt_columns = ["rt_t_values[$i]" for i in 0:(length(alpha_times)-1)]
other_columns = [:gamma, :nu, :e0, :i0, :rw_sigma]
# rewrite other_columns as characters
other_columns = [string(col) for col in other_columns]
all_columns = vcat(other_columns, rt_columns, ["log_likelihood", "actual_iteration"])
my_samples_frame = DataFrame(my_samples, all_columns)
# convert all_columns into symbols
all_columns = [Symbol(col) for col in all_columns]
chain = Chains(my_samples, all_columns)
# Calculate the effective sample size
my_ess = ess(chain)
sim = "LBR_version2"
seed = 123467
CSV.write(resultsdir("my_generated_quantities", 
string("ess_joint_hetchron_inseq_samples_sim", sim, "_seed", seed, ".csv")), my_samples_frame)
