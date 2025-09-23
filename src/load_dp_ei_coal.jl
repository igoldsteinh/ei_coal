# file for loading priors and data for the ei coal model 
if sim_id == 1
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(2.0)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = []
end
if sim_id == 2
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(2.0)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = []
end
if sim_id == 3
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(2.0)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = ones(Int, 20)
end
if sim_id == 4
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(2.0)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = ones(Int, 30)
end
if sim_id == 5
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(1.2)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = []
end
if sim_id == 6
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(1.2)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = []
end
if sim_id == 7
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(1)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = ones(Int, 20)
end
if sim_id == 8
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(1.2)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = ones(Int, 30)
end
if sim_id == 9
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(2.0)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = []
end
if sim_id == 10
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(2.0)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = []
end
if sim_id == 11
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(2.0)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = ones(Int, 20)
end
if sim_id == 12
  sim_dict = CSV.read("data/sim_data/sim_dict.csv", DataFrame)
  # filter sim_dict by sim 
  sim_vals = filter(row -> row[:sim_id] == sim_id, sim_dict)
  all_trees = CSV.read(datadir("sim_data", sim_vals.tree_file[1]), DataFrame)
  ## Define Priors
  const log_rt_init_mean = log(2.0)
  const log_rt_init_sd = 0.2
  const log_gamma_mean = log(1/4)
  const log_gamma_sd = 0.25
  const log_nu_mean = log(1/7)
  const log_nu_sd = 0.25
  const log_e0_mean = log(1.1)
  const log_e0_sd = 0.05
  const log_i0_mean = log(1.1)
  const log_i0_sd = 0.05
  const log_rw_mean = log(0.2)
  const log_rw_sigma_sd = 0.1
  log_prior_means = [log_gamma_mean, log_nu_mean, log_e0_mean, log_i0_mean,  log_rw_mean, log_rt_init_mean,]
  cholC = create_cholc_matrix(log_gamma_sd, log_nu_sd, log_e0_sd, log_i0_sd, log_rw_sigma_sd, log_rt_init_sd)
  tree = DataFrames.subset(all_trees, :sim => ByRow(x -> x == sim_num))
  reverse_samp_lin = ones(Int, 30)
end
