# PARTICLES-IN-A-POTENTIAL

This repository contains the code that was used to simulate and analyse a model of a many particle system. The goal was to sample with Markov Chain Monte Carlo from the density of states and use those sample to find the phases and phase transitions of the system. 

The density of states is given by rho = 1 / Z * exp[-1/T * E(X)] with Z the partition function (unimportant normalization constant in this case) and the total energy E(X) of the system given by the pair potential V_{ij} as E(X) = sum_{i<j} V_{ij} (r_i, r_j). The pair potential is assumed to be V_{ij} = q_i * q_j / |r_i - r_j| + 1 / |r_i - r_j|^8.

The code is structured as follows: Simulation and analysis is done inside the R Markdown steering file steering.Rmd. The steering file loads the R scripts DataFrameHelper.R, MCMC_functions.R and MCMC_functions.cpp (contains Rcpp code to use in the R environment). The loaded R scripts provide all necessary functions to generate samples of the system.
