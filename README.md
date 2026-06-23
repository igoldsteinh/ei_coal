# ei_coal
This is the code repository for the paper Coalescent Inference for Epidemics with Latent Periods.
The main model was fit in `Julia`, while visualization of results was done in `R`. 
## Navigation
```
├── data                          <- Processed real and simulated data
│   ├── compare_data              <- Data for comparing coalescent simulators
│   ├── sim_data                  <- Simulated data for model fitting
|   └── real_data                 <- data used for Liberia Ebola analysis
|    
│
├── figures                       <- Paper figures
│
├── results                       <- Model outputs
│   ├── my_generated_quantities   <- Correctly scaled posteriors, quantiles and mcmc samples
│   ├── my_mcmc_summaries         <- mcmc diagnostic summaries
│
├── scripts                       <- Paper code 
│   ├── BEAST2                    <- BDMM xml, processing files, results
│   ├── ebola_scripts             <- Process real data, BEAST X .xml file, fit model julia file
│   ├── phydynR                   <- fit PhydynR, process results, results
│   ├── process_results           <- Final processing of mcmc/summarising simulation results
│   ├── simulated_data            <- Simulate data
│   └── visualize_results         <- Turn summaries of model results into paper figures
│   
├── src                           <- Models, priors, simulation engines, utility functions
```

## Setting up the `Julia` environment. 
The results from this project were generated using `Julia 1.11.4`. 
We recommend installing `Julia` via [`juliaup`](https://github.com/JuliaLang/juliaup).
Once you have `Julia` installed, from the terminal, navigate to the project root directory then type `julia`. 
Your terminal will look like:
```
julia>
```
Now type `]`. Your terminal should now look like:
```
(@v1.11.4) pkg>
```
Then use the following commands
```
activate .
```
and 
```
instantiate
```
More information on `Julia` environments is available in the [Environments documentation](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project).

## Simulation name key
When executing scripts, the `sim_id` parameter controls what simulation is being used, the `sim_num` parameter controls the seed and also the specific data set used. 
For simulations, we used values of `sim_num` from 1 to 100. 
Here is a key translating the values of `sim_id`:
* `sim_id=1` = `Fixed, N=50, Isochronous`
* `sim_id=2` = `Fixed, N=100, Isochronous`
* `sim_id=3` = `Fixed, N=50, Heterochronous`
* `sim_id=4` = `Fixed, N=100, Heterochronous`
* `sim_id=5` = `Increase, N=50, Isochornous`
* `sim_id=6` = `Increase, N=100, Isochronous`
* `sim_id=7` = `Increase, N=50, Heterechronous`
* `sim_id=8` = `Increase, N=100, Heterochronous`
* `sim_id=9` = `Control, N=50, Isochronous`
* `sim_id=10` = `Control, N=100, Isochronous`
* `sim_id=11` = `Control, N=50, Heterochronous`
* `sim_id=12` = `Control, N=100, Heterochronous`
More detailed information on simulations are stored in the [simulation key file](https://github.com/igoldsteinh/ei_coal/blob/main/data/sim_data/sim_dict.csv).
