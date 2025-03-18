# Metaautomaton
Scripts and data for reproducing results from "Assortative dispersal facilitates the regional maintenance of alternative stable states" 

## `functions`
- all functions required for simulating the model are within the `final_funs.R` script

## `simulation scripts`
- Small landscape simulations are split into two: $\alpha_{ij} > 1$ (`ass-sims`) and $\alpha_{ij} \leq 1$ (`coex-sims`) 
- Simulations for the landscape portraits (`landscape-portaits`) include only information about the state space (has no other metrics)
- Simulations to generate the stationary distribution of metrics are contained within `stationary-states`
- Simulations for the time-series of macroscopic variables are contained within `landscape-timeseries`
- All simulations are ran in parallel on HPCs with at least 25 cores and at most 120 cores 
