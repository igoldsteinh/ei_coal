# ei_coal
This is the code repository for the paper Coalescent Inference for Epidemics with Latent Periods.
The main model was fit in `Julia`, while visualization of results was done in `R`. 
## Navigation
```
├── data                          <- Processed real and simulated data
│   ├── sim_data                  <- Simulated data for model fitting
|   └── compare_data              <- Data for comparing coalescent simulators 
│
├── figures                       <- Paper figures
│
├── results                       <- Model outputs
│   ├── my_generated_quantities   <- Correctly scaled posteriors, quantiles and mcmc samples
│   ├── my_mcmc_summaries         <- mcmc diagnostic summaries
│
├── scripts                       <- Paper code 
│   ├── ebola_scripts             <- Process real data, BEAST X .xml file, fit model julia file
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
When executing scripts, the `sim` parameter controls what simulation is being used, the `seed` parameter controls the seed and also the specific data set used. 
For simulations, we used values of `seed` from 1 to 100. 
Here is a key translating the values of `sim`:
* `sim=1` = `Fixed, N=50, Isochronous`
* `sim=2` = `Fixed, N=100, Isochronous`
* `sim=3` = `Fixed, N=50, Heterochronous`
* `sim=4` = `Fixed, N=100, Heterochronous`
* `sim=5` = `Increase, N=50, Isochornous`
* `sim=6` = `Increase, N=100, Isochronous`
* `sim=7` = `Increase, N=50, Heterechronous`
* `sim=8` = `Increase, N=100, Heterochronous`
* `sim=9` = `Control, N=50, Isochronous`
* `sim=10` = `Control, N=100, Isochronous`
* `sim=11` = `Control, N=50, Heterochronous`
* `sim=12` = `Control, N=100, Heterochronous`
More detailed information on simulations are stored in the [simulation key file](https://github.com/igoldsteinh/ei_coal/blob/main/data/sim_data/sim_dict.csv).