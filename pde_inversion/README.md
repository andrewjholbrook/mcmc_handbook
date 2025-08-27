# Bayesian Inversion of Nonlinear Partial Differential Equations
This directory contains code for the "Bayesian Inversion of Nonlinear Partial Differential Equations" example in the Glatt-Holtz/Holbrook/Krometis/Mondaini/Sheth chapter [1] of the 2025 version of the _Handbook of Markov Chain Monte Carlo_.

It relies on three other repositories of mine:
- MCMC implementations from here: https://github.com/krometis/InfDimMCMC.jl
- Advection-Diffusion equation solvers from here: https://github.com/krometis/AdvectionDiffusion.jl
- Spectral discretization used by the solvers above from here: https://github.com/krometis/SpectralDiscrete2D.jl

Chains can be run by calling `run.jl` for a given scenario. Traditional MCMC methods (e.g., pCN or HMC) can be run using the scenario associated with the original 2020 JUQ paper [2], e.g.,
```
julia scenarios/2020juqex2/run.jl --scen=2020juqex2 --nsamp=250000 --mcmc="pcn|0.2" --ar=0.0
```

Neural network-based samplers can be run using the `lolmcNN`* scenario, e.g.,
```
julia scenarios/lolmcNN/run.jl --scen=lolmcNN --nsamp=100000 --mcmc="lolhmc|0.25|2" --ar=0.0 --modelfile="trained_models/large.bson"
```

*`lol` stands for "learning-optimized Langevin" and was an old joke amongst ourselves - sorry you're now being subjected to it. :-/

Trained neural networks have been provided in the `trained_models` section; they were trained using the script `scenarios/lolmcNN/trainNN.jl`.


## References
[1] Glatt-Holtz, N. E.; Holbrook, A. J.; Krometis, J. A.; Mondaini, C. F.; Sheth, A. Sacred and Profane: From the Involutive Theory of MCMC to Helpful Hamiltonian Hacks. arXiv October 22, 2024. https://doi.org/10.48550/arXiv.2410.17398.

[2] Borggaard, J.; Glatt-Holtz, N.; Krometis, J. A Bayesian Approach to Estimating Background Flows from a Passive Scalar. SIAM/ASA J. Uncertainty Quantification 2020, 8 (3), 1036–1060. https://doi.org/10.1137/19M1267544.
