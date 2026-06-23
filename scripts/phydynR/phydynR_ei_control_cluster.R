# Fit phydynR to simulated data
# this script is based on the tutorial http://emvolz-phylodynamics.github.io/phydynR/articles/SenegalHIVmodel.html
# by Fabricia F Nascimiento
library(ape)
library(akima)
library(BayesianTools)
library(coda)
library(phydynR)
library(treedater)
library(ggplot2)
library(readr)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  sim_num_val = 6
} else {
  sim_num_val <- as.integer(args[1])
}

# define our deme matrix
ei_demes <- c("E", "I")
ei_births <- matrix(nrow = 2, ncol = 2)
rownames(ei_births) <- ei_demes
colnames(ei_births) <- ei_demes

ei_births["I", "E"] <- "parms$rt_traj(t, parms) * parms$nu * I"
ei_births["E", "I"] <- "0."
ei_births["E", "E"] <- "0."
ei_births["I", "I"] <- "0."

# migration matrix
ei_migs <- matrix(nrow = 2, ncol = 2)
rownames(ei_migs) <- ei_demes
colnames(ei_migs) <- ei_demes
ei_migs["E", "I"] <- "parms$gamma * E"
ei_migs["I", "E"] <- "0."
ei_migs["E", "E"] <- "0."
ei_migs["I", "I"] <- "0."

# death matrix
ei_deaths <- c("0.", "parms$nu * I")
names(ei_deaths) <- ei_demes

# define parameters
alpha_times <- c(0, seq(from = 7.0, by = 7.0, length.out = 21))


ei_THETA <- list(
  inite0 = 1.1,
  initi0 = 1.1,
  gamma = 1/4,
  nu = 1/7,
  rt_traj = function(t, parms) {
    approx(
      x = alpha_times,
      y = parms$rt_values,
      xout = t,
      method = "constant",
      f = 0,  
      rule = 2
    )$y
  },
  rt_values = rep(2, 22),
  rw_sd = 0.2
)

ei_dm <- build.demographic.process(births = ei_births,
                                deaths = ei_deaths,
                                migrations = ei_migs,
                                parameterNames = names(ei_THETA),
                                rcpp = FALSE,
                                sde = FALSE)

# load tree data
ei_tree <- read.tree(here::here("data", 
                                 "sim_data", 
                                 "control_trees",
                                 paste0("control_iso50", "_simnum", sim_num_val, ".tree")))

# meta data
E <- rep(0, length(ei_tree$tip.label))
I <- rep(1,  length(ei_tree$tip.label))

ei_sampleStates <- cbind(E, I)
rownames(ei_sampleStates) <- ei_tree$tip.label

# sample times
tree_height <- max(node.depth.edgelength(ei_tree))
min_length <- min(ei_tree[["edge.length"]])

ei_times <- rep(ceiling(tree_height), length(ei_tree$tip.label))
names(ei_times) <- ei_tree$tip.label
# create dated tree object in phylodyn
ei_dated.tree <- phydynR::DatedTree(phylo = ei_tree,
                                 sampleTimes = ei_times,
                                 sampleStates = ei_sampleStates,
                                 minEdgeLength = min_length,
                                 tol = 0.1)

# likelihood calculation
# set initial compartment counts
ei_X0 = c(E = 1.1,
          I = 1.1)
phydynR::colik(tree = ei_dated.tree,
               theta = ei_THETA,
               demographic.process.model = ei_dm,
               x0 = ei_X0,
               t0 = 0,
               res = 1e3,
               timeOfOriginBoundaryCondition = FALSE,
               AgtY_penalty = 1,
               maxHeight = 160)

# parameter inference object creation
obj_fun <- function(parameters){
  parameters <- unname(parameters)
  
  # add the values of THETA to a new variable named THETA.new
  THETA.new <- ei_THETA
  
  # change the values in THETA.new to the new proposals that will be evaluated
  THETA.new$inite0 <- parameters[1]
  THETA.new$initi0 <- parameters[2]
  THETA.new$gamma <- parameters[3]
  THETA.new$nu <- parameters[4]
  THETA.new$rt_values <- parameters[5:26]
  THETA.new$rw_sd <- parameters[27]

  X0 <- c(E = unname(THETA.new$inite0),
          I = unname(THETA.new$initi0))
  
  # After changing the parameter values to the new proposals, a likelihood is
  # calculated with the function phydynR::colik.
  mll <- colik(tree = ei_dated.tree,
               theta = THETA.new,
               demographic.process.model = ei_dm,
               x0 = X0,
               t0 = 0,
               res = 1e3, #TODO
               timeOfOriginBoundaryCondition = FALSE,
               AgtY_penalty = 1,
               maxHeight = 160)
  
  return(mll)
  
}

# prior specification
# phydynR has no default random walk prior
# but it can be created using the conditional form of the random walk
densities <-  function(par){
  if (any(par < 0)) {return -Inf}
  d1 = dlnorm(par[1], meanlog = log(1.1), sdlog = 0.05, log = TRUE) #inite0
  d2 = dlnorm(par[2], meanlog = log(1.1), sdlog = 0.05, log = TRUE) #initi0
  d3 = dlnorm(par[3], meanlog = log(1/4), sdlog = 0.2, log = TRUE) #gamma
  d4 = dlnorm(par[4], meanlog = log(1/7), sdlog = 0.2, log = TRUE) #nu
  d5 = dlnorm(par[5], meanlog = log(2), sdlog = 0.2, log = TRUE) #init rt
  # next the rest of the random walk 
  d6 = dlnorm(par[6], meanlog = log(par[5]), sdlog = par[27], log = TRUE)
  d7 = dlnorm(par[7], meanlog = log(par[6]), sdlog = par[27], log = TRUE)
  d8 = dlnorm(par[8], meanlog = log(par[7]), sdlog = par[27], log = TRUE)
  d9 = dlnorm(par[9], meanlog = log(par[8]), sdlog = par[27], log = TRUE)
  d10 = dlnorm(par[10], meanlog = log(par[9]), sdlog = par[27], log = TRUE)
  d11 = dlnorm(par[11], meanlog = log(par[10]), sdlog = par[27], log = TRUE)
  d12 = dlnorm(par[12], meanlog = log(par[11]), sdlog = par[27], log = TRUE)
  d13 = dlnorm(par[13], meanlog = log(par[12]), sdlog = par[27], log = TRUE)
  d14 = dlnorm(par[14], meanlog = log(par[13]), sdlog = par[27], log = TRUE)
  d15 = dlnorm(par[15], meanlog = log(par[14]), sdlog = par[27], log = TRUE)
  d16 = dlnorm(par[16], meanlog = log(par[15]), sdlog = par[27], log = TRUE)
  d17 = dlnorm(par[17], meanlog = log(par[16]), sdlog = par[27], log = TRUE)
  d18 = dlnorm(par[18], meanlog = log(par[17]), sdlog = par[27], log = TRUE)
  d19 = dlnorm(par[19], meanlog = log(par[18]), sdlog = par[27], log = TRUE)
  d20 = dlnorm(par[20], meanlog = log(par[19]), sdlog = par[27], log = TRUE)
  d21 = dlnorm(par[21], meanlog = log(par[20]), sdlog = par[27], log = TRUE)
  d22 = dlnorm(par[22], meanlog = log(par[21]), sdlog = par[27], log = TRUE)
  d23 = dlnorm(par[23], meanlog = log(par[22]), sdlog = par[27], log = TRUE)
  d24 = dlnorm(par[24], meanlog = log(par[23]), sdlog = par[27], log = TRUE)
  d25 = dlnorm(par[25], meanlog = log(par[24]), sdlog = par[27], log = TRUE)
  d26 = dlnorm(par[26], meanlog = log(par[25]), sdlog = par[27], log = TRUE)
  d27 = dlnorm(par[27], meanlog = log(0.2), sdlog = 0.1, log = TRUE)
  return(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9 + d10 + d11 + d12 + d13 + d14 + d15 +
         d16 + d17 + d18 + d19 + d20 + d21 + d22 + d23 + d24 + d25 + d26 + d27)
}

# bayesian tools sampler
sampler <-  function(n=1){
  d1 = rlnorm(n, meanlog = log(1.1), sdlog = 0.05) #inite0
  d2 = rlnorm(n, meanlog = log(1.1), sdlog = 0.05) #initi0
  d3 = rlnorm(n, meanlog = log(1/4), sdlog = 0.2) #gamma
  d4 = rlnorm(n, meanlog = log(1/7), sdlog = 0.2) #nu
  d5 = rlnorm(n, meanlog = log(2), sdlog = 0.2) #init rt
  d27 = rlnorm(n, meanlog = log(0.2), sdlog = 0.1) #rw prior
  other_rts = rep(0.0, length = 21)
  other_rts[1] = rlnorm(n, meanlog = log(d5), sdlog = d27)
  for (i in 2:21) {
    other_rts[i] = rlnorm(n, meanlog = log(other_rts[i-1]), sdlog = d27)
    
  }
  return(cbind(d1, d2, d3, d4, d5, t(other_rts), d27))
}
# Create prior (necessary for the BayesianTools package)
prior <- createPrior(density = densities,
                     sampler = sampler)

# more bayesian set up objects

settings = list(iterations = 10000, nrChains = 1, thin = 1)
# Create bayesianSetup
bayesianSetup <- createBayesianSetup(likelihood = obj_fun , prior = prior)
# run the MCMC
start_time_pre = Sys.time()
out <- BayesianTools::runMCMC(bayesianSetup = bayesianSetup,
                              sampler = "DEzs",
                              settings = settings)
end_time_pre = Sys.time()
print(end_time_pre - start_time_pre)

out_sample <- BayesianTools::getSample(out, start = 1000, coda = TRUE)

# continuing with new z matrix
x <- BayesianTools::getSample(out, start = 20)

# Get the range for the parameter estimates for the previous run
rangePost = apply(x, 2, range)

#get unique values of x
u_x <- unique(x)

#create new Z matrix based on previous run
num_vals = 120 * 27
newZ = matrix(runif(num_vals, rangePost[1,], rangePost[2,]), ncol = 27, byrow = T)

# run chains in parallel
pos1 = 72
pos2 = 73
pos3 = 74
pos4 = 80
iter = 80000 # number of iterations
settings = list(Z = newZ, 
                startValue =  u_x[c(pos1, pos2, pos3, pos4), ], 
                nrChains = 1, 
                iterations = iter, 
                thin = 10)


# Create bayesianSetup
bayesianSetup <- createBayesianSetup(likelihood = obj_fun,
                                    prior = prior,
                                    parallel = 4)

set.seed(sim_num_val)
start_time = Sys.time()
outZ <- runMCMC(bayesianSetup = bayesianSetup,  sampler = "DEzs", settings = settings )
end_time = Sys.time()
stopParallel(bayesianSetup)
print(end_time - start_time)

out_sample <- BayesianTools::getSample(outZ, start = 1, coda = TRUE)

# create outputs
mcmc_matrix_coda <- as.matrix(out_sample, chains = TRUE)
mcmc_df <- as.data.frame(mcmc_matrix_coda)
new_names <- c(".chain", "e0", "i0", "gamma", "nu", "rt_t_values[0]",  "rt_t_values[1]",  "rt_t_values[2]",
               "rt_t_values[3]", "rt_t_values[4]",  "rt_t_values[5]",  "rt_t_values[6]",  
               "rt_t_values[7]", "rt_t_values[8]",  "rt_t_values[9]", 
               "rt_t_values[10]", "rt_t_values[11]", "rt_t_values[12]", 
               "rt_t_values[13]", "rt_t_values[14]", "rt_t_values[15]", "rt_t_values[16]",
               "rt_t_values[17]", "rt_t_values[18]", "rt_t_values[19]", "rt_t_values[20]",
               "rt_t_values[21]", "rw_sigma")
names(mcmc_df) <- new_names
write_csv(mcmc_df, here::here("scripts", "phydynR", 
                              paste0("phydynR_controliso50_simnum", sim_num_val, "_eires.csv")))