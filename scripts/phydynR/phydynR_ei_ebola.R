# phydynR fit to ebola liberia 
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
set.seed(1235689)
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
alpha_times <- c(0, seq(from = 7.0, by = 7.0, length.out = 52))


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
  rt_values = rep(2, 53),
  rw_sd = 0.2
)

ei_dm <- build.demographic.process(births = ei_births,
                                deaths = ei_deaths,
                                migrations = ei_migs,
                                parameterNames = names(ei_THETA),
                                rcpp = FALSE,
                                sde = FALSE)

# load tree data
ei_tree <- read.nexus(here::here("data", "real_data", "liberia_mcctreev2.nexus"))
ei_tree[["edge.length"]] <- ei_tree[["edge.length"]] * 365.25
# meta data
E <- rep(0, length(ei_tree$tip.label))
I <- rep(1,  length(ei_tree$tip.label))
ei_sampleStates <- cbind(E, I)
rownames(ei_sampleStates) <- ei_tree$tip.label

# Get the distance from root to each tip
tip_depths <- node.depth.edgelength(ei_tree)
# Extract only the tips (first n elements, where n = number of tips)
n_tips <- length(ei_tree$tip.label)
ei_times <- tip_depths[1:n_tips]
# Name the vector with tip labels
names(ei_times) <- ei_tree$tip.label
tree_height <- max(node.depth.edgelength(ei_tree))
min_length <- min(ei_tree[["edge.length"]])
# created dated tree object in phylodyn
ei_dated.tree <- phydynR::DatedTree(phylo = ei_tree,
                                 sampleTimes = ei_times,
                                 sampleStates = ei_sampleStates,
                                 minEdgeLength = min_length,
                                 tol = 0.1)

# likelihood calculation
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
               maxHeight = 366)

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
  THETA.new$rt_values <- parameters[5:57]
  THETA.new$rw_sd <- parameters[58]
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
               maxHeight = 366)
  
  return(mll)
  
}

# prior specification
densities <-  function(par){
  if (any(par < 0)) {return -Inf}
  d1 = dlnorm(par[1], meanlog = log(1.1), sdlog = 0.05, log = TRUE) #inite0
  d2 = dlnorm(par[2], meanlog = log(1.1), sdlog = 0.05, log = TRUE) #initi0
  d3 = dlnorm(par[3], meanlog = log(1/7), sdlog = 0.45, log = TRUE) #gamma
  d4 = dlnorm(par[4], meanlog = log(1/7), sdlog = 0.3, log = TRUE) #nu
  d5 = dlnorm(par[5], meanlog = 0.7, sdlog = 0.5, log = TRUE) #init rt
  # next the rest of the random walk 
  d6 = dlnorm(par[6], meanlog = log(par[5]), sdlog = par[58], log = TRUE)
  d7 = dlnorm(par[7], meanlog = log(par[6]), sdlog = par[58], log = TRUE)
  d8 = dlnorm(par[8], meanlog = log(par[7]), sdlog = par[58], log = TRUE)
  d9 = dlnorm(par[9], meanlog = log(par[8]), sdlog = par[58], log = TRUE)
  d10 = dlnorm(par[10], meanlog = log(par[9]), sdlog = par[58], log = TRUE)
  d11 = dlnorm(par[11], meanlog = log(par[10]), sdlog = par[58], log = TRUE)
  d12 = dlnorm(par[12], meanlog = log(par[11]), sdlog = par[58], log = TRUE)
  d13 = dlnorm(par[13], meanlog = log(par[12]), sdlog = par[58], log = TRUE)
  d14 = dlnorm(par[14], meanlog = log(par[13]), sdlog = par[58], log = TRUE)
  d15 = dlnorm(par[15], meanlog = log(par[14]), sdlog = par[58], log = TRUE)
  d16 = dlnorm(par[16], meanlog = log(par[15]), sdlog = par[58], log = TRUE)
  d17 = dlnorm(par[17], meanlog = log(par[16]), sdlog = par[58], log = TRUE)
  d18 = dlnorm(par[18], meanlog = log(par[17]), sdlog = par[58], log = TRUE)
  d19 = dlnorm(par[19], meanlog = log(par[18]), sdlog = par[58], log = TRUE)
  d20 = dlnorm(par[20], meanlog = log(par[19]), sdlog = par[58], log = TRUE)
  d21 = dlnorm(par[21], meanlog = log(par[20]), sdlog = par[58], log = TRUE)
  d22 = dlnorm(par[22], meanlog = log(par[21]), sdlog = par[58], log = TRUE)
  d23 = dlnorm(par[23], meanlog = log(par[22]), sdlog = par[58], log = TRUE)
  d24 = dlnorm(par[24], meanlog = log(par[23]), sdlog = par[58], log = TRUE)
  d25 = dlnorm(par[25], meanlog = log(par[24]), sdlog = par[58], log = TRUE)
  d26 = dlnorm(par[26], meanlog = log(par[25]), sdlog = par[58], log = TRUE)
  d27 = dlnorm(par[27], meanlog = log(par[26]), sdlog = par[58], log = TRUE)
  d28 = dlnorm(par[28], meanlog = log(par[27]), sdlog = par[58], log = TRUE)
  d29 = dlnorm(par[29], meanlog = log(par[28]), sdlog = par[58], log = TRUE)
  d30 = dlnorm(par[30], meanlog = log(par[29]), sdlog = par[58], log = TRUE)
  d31 = dlnorm(par[31], meanlog = log(par[30]), sdlog = par[58], log = TRUE)
  d32 = dlnorm(par[32], meanlog = log(par[31]), sdlog = par[58], log = TRUE)
  d33 = dlnorm(par[33], meanlog = log(par[32]), sdlog = par[58], log = TRUE)
  d34 = dlnorm(par[34], meanlog = log(par[33]), sdlog = par[58], log = TRUE)
  d35 = dlnorm(par[35], meanlog = log(par[34]), sdlog = par[58], log = TRUE)
  d36 = dlnorm(par[36], meanlog = log(par[35]), sdlog = par[58], log = TRUE)
  d37 = dlnorm(par[37], meanlog = log(par[36]), sdlog = par[58], log = TRUE)
  d38 = dlnorm(par[38], meanlog = log(par[37]), sdlog = par[58], log = TRUE)
  d39 = dlnorm(par[39], meanlog = log(par[38]), sdlog = par[58], log = TRUE)
  d40 = dlnorm(par[40], meanlog = log(par[39]), sdlog = par[58], log = TRUE)
  d41 = dlnorm(par[41], meanlog = log(par[40]), sdlog = par[58], log = TRUE)
  d42 = dlnorm(par[42], meanlog = log(par[41]), sdlog = par[58], log = TRUE)
  d43 = dlnorm(par[43], meanlog = log(par[42]), sdlog = par[58], log = TRUE)
  d44 = dlnorm(par[44], meanlog = log(par[43]), sdlog = par[58], log = TRUE)
  d45 = dlnorm(par[45], meanlog = log(par[44]), sdlog = par[58], log = TRUE)
  d46 = dlnorm(par[46], meanlog = log(par[45]), sdlog = par[58], log = TRUE)
  d47 = dlnorm(par[47], meanlog = log(par[46]), sdlog = par[58], log = TRUE)
  d48 = dlnorm(par[48], meanlog = log(par[47]), sdlog = par[58], log = TRUE)
  d49 = dlnorm(par[49], meanlog = log(par[48]), sdlog = par[58], log = TRUE)
  d50 = dlnorm(par[50], meanlog = log(par[49]), sdlog = par[58], log = TRUE)
  d51 = dlnorm(par[51], meanlog = log(par[50]), sdlog = par[58], log = TRUE)
  d52 = dlnorm(par[52], meanlog = log(par[51]), sdlog = par[58], log = TRUE)
  d53 = dlnorm(par[53], meanlog = log(par[52]), sdlog = par[58], log = TRUE)
  d54 = dlnorm(par[54], meanlog = log(par[53]), sdlog = par[58], log = TRUE)
  d55 = dlnorm(par[55], meanlog = log(par[54]), sdlog = par[58], log = TRUE)
  d56 = dlnorm(par[56], meanlog = log(par[55]), sdlog = par[58], log = TRUE)
  d57 = dlnorm(par[57], meanlog = log(par[56]), sdlog = par[58], log = TRUE)
  d58 = dlnorm(par[58], meanlog = log(0.05), sdlog = 0.2, log = TRUE)
  return(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9 + d10 + d11 + d12 + d13 + d14 + d15 +
           d16 + d17 + d18 + d19 + d20 + d21 + d22 + d23 + d24 + d25 + d26 + d27 + 
           d28 + d29 + d30 + d31 + d32 + d33 + d34 + d35 + d36 + d37 + d38 + d39 + 
           d40 + d41 + d42 + d43 + d44 + d45 + d46 + d47 + d48 + d49 + d50 + d51 + 
           d52 + d53 + d54 + d55 + d56 + d57 + d58)
}

# bayesian tools sampler
sampler <-  function(n=1){
  d1 = rlnorm(n, meanlog = log(1.1), sdlog = 0.05) #inite0
  d2 = rlnorm(n, meanlog = log(1.1), sdlog = 0.05) #initi0
  d3 = rlnorm(n, meanlog = log(1/7), sdlog = 0.45) #gamma
  d4 = rlnorm(n, meanlog = log(1/7), sdlog = 0.3) #nu
  d5 = rlnorm(n, meanlog = 0.7, sdlog = 0.5) #init rt
  d58 = rlnorm(n, meanlog = log(0.05), sdlog = 0.2) #rw prior
  other_rts = rep(0.0, length = 52)
  other_rts[1] = rlnorm(n, meanlog = log(d5), sdlog = d58)
  for (i in 2:52) {
    other_rts[i] = rlnorm(n, meanlog = log(other_rts[i-1]), sdlog = d58)
    
  }
  return(cbind(d1, d2, d3, d4, d5, t(other_rts), d58))
}
test = sampler(n=1)
densities(test)
# create the prior
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
#cretae new Z matrix based on previous run
num_vals = 120 * 58
newZ = matrix(runif(num_vals, rangePost[1,], rangePost[2,]), ncol = 58, byrow = T)
# run iterations in parallel
pos1 = 13
pos2 = 23
pos3 = 45
pos4 = 76
iter = 200000 # number of iterations
settings = list(Z = newZ, 
                startValue =  u_x[c(pos1, pos2, pos3, pos4), ], 
                nrChains = 1, 
                iterations = iter, 
                thin = 1)

# Create bayesianSetup
bayesianSetup <- createBayesianSetup(likelihood = obj_fun,
                                    prior = prior,
                                    parallel = 4)

set.seed(1234)
start_time = Sys.time()
outZ <- runMCMC(bayesianSetup = bayesianSetup,  sampler = "DEzs", settings = settings )
end_time = Sys.time()
stopParallel(bayesianSetup)
print(end_time - start_time)
library(readr)
out_sample <- BayesianTools::getSample(outZ, start = 1, coda = TRUE)
# process results
mcmc_matrix_coda <- as.matrix(out_sample, chains = TRUE)
mcmc_df <- as.data.frame(mcmc_matrix_coda)
new_names <- c(".chain", "e0", "i0", "gamma", "nu", 
               "rt_t_values[0]",  "rt_t_values[1]",  "rt_t_values[2]",
               "rt_t_values[3]",  "rt_t_values[4]",  "rt_t_values[5]",  
               "rt_t_values[6]",  "rt_t_values[7]",  "rt_t_values[8]",  
               "rt_t_values[9]",  "rt_t_values[10]", "rt_t_values[11]", 
               "rt_t_values[12]", "rt_t_values[13]", "rt_t_values[14]", 
               "rt_t_values[15]", "rt_t_values[16]", "rt_t_values[17]", 
               "rt_t_values[18]", "rt_t_values[19]", "rt_t_values[20]",
               "rt_t_values[21]", "rt_t_values[22]", "rt_t_values[23]",
               "rt_t_values[24]", "rt_t_values[25]", "rt_t_values[26]",
               "rt_t_values[27]", "rt_t_values[28]", "rt_t_values[29]",
               "rt_t_values[30]", "rt_t_values[31]", "rt_t_values[32]",
               "rt_t_values[33]", "rt_t_values[34]", "rt_t_values[35]",
               "rt_t_values[36]", "rt_t_values[37]", "rt_t_values[38]",
               "rt_t_values[39]", "rt_t_values[40]", "rt_t_values[41]",
               "rt_t_values[42]", "rt_t_values[43]", "rt_t_values[44]",
               "rt_t_values[45]", "rt_t_values[46]", "rt_t_values[47]",
               "rt_t_values[48]", "rt_t_values[49]", "rt_t_values[50]",
               "rt_t_values[51]", "rt_t_values[52]",
               "rw_sigma")
names(mcmc_df) <- new_names
write_csv(mcmc_df, here::here("scripts", "phydynR", "phydynR_ebola_res.csv"))
# note this file is very large, so it is not stored on the repo