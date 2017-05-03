# Sample size calculation for Stratified Micro-Randomized Trial

The below table provides and overview of the R scripts.
The R scripts are split by case (conditional or marginal) and then again by type (overall or robustness testing).
R output and simulation data are respectively saved to .Rout and .RData files by the same name.

File | Description
---- | ----
[functions.R](functions.R) | All functions related to simulation and estimation are stored in this file
[setup.R](setup.R) | File which contains all the necessary setup details
[estimation.R](estimation.R) | File which runs the simulation study (in parallel -- setup for Flux use) given setup and sample size
