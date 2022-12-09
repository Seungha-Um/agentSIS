# Particle MCMC for agent-based SIS model

This repository contains code to implement the methodology described in
the paper “Bayesian Inference on Agent-based Models of Disease
Transmission”, by Um and Adhikari (2022+)

This package uses the primary functions from
[`agents`](https://github.com/nianqiaoju/agents).

# Installation

The packages can be installed with the `devtools` package:

    library(devtools) 
    devtools::install_github(repo='Seungha-Um/agentSIS') 

The package with the vignettes can be installed with

    devtools::install_github(repo='Seungha-Um/agentSIS', build_vignettes = TRUE) 

and then accessed by running `browseVignettes("agentSIS")` (to reproduce
our results, one will need to increase the number of MCMC samples).
Alternatively, vignettes are available at
[Simulation](https://rpubs.com/sheom0808/981711)
