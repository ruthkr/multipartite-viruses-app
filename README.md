# Multipartite Viruses

These scripts were created as a part of my master thesis (with [Tomas Alarcon](https://www.icrea.cat/Web/ScientificStaff/tomas-alarcon-210760) as my supervisor and [Josep Sardanyes](http://complex.upf.es/~josep/Site/Welcome.html) as my co-supervisor. The scripts are used to analysis the determinsitic models for the behaviour of multipartite viruses. More info about multipartie viruses : [A. Sicard, et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5094692/) and [A. Lucía Sanz, et al](https://www.nature.com/articles/s41540-017-0035-y).

There are two main parts of this project, which are deterministic model and stochastic model for each virus replication model and full cycle bipartite viruses. 

## Replication Model
### Deterministic analysis
- `2eq_rk4_bipartite_test.c`: simulate the replication model using forth order Runge-Kutta
- `2eq_rk4_bipartite_var_gamma.c`: bifurcation analysis as the paramater **gamma** changes
- `2eq_rk4_bipartite_var_kappa.c`: bifurcation analysis as the paramater **kappa** changes
- `2eq_rk4_rep_var_gamma_tozero.c`: estimate the separatrix surface of the model
- `2eq_plot.R`: plot the simulation results from `2eq_rk4_bipartite_test.c`
- `2eq_phase_portrait.R`: plot of simple phase portrait
- `2eq_plot_gamma_var.R`: plot the bifurcation diagram for parameter **gamma**
- `2eq_plot_kappa_var.R`: plot the bifurcation diagram for parameter **kappa**


### Stochastic analysis
- `2eq_ssa.c`: simulate the Markov chain model of replication model using Gillespie's algorithm (python version `2eq_ssa.py`)
- `2eq_ssa_prob_exit.c`: calculate the probability of extinction due to intrinsic noise
- `2eq_plot_gillespie.R`: plot the simulation results of CTCM model
- `2eq_separatrix_point.R`: plot deteministic and stochastic trajectories

## Bipartite Model
### Deterministic analysis
- `6eq_rk4_bipartite_test.c`: simulate the bipartite viruses' model using forth order Runge-Kutta
- `6eq_rk4_bipartite_var_gamma.c`: simulate the bifurcation analysis for parameter **gamma**
- `6eq_plot.R`:plot the simulation results from `6eq_rk4_bipartite_test.c`
- `6eq_ratio_plot_gamma.R`: plot the ratio of the simulation going to origin over all of the total simulation
- `6eq_plot_gamma_var.R`: plot the bifurcation diagram for parameter **gamma**

### Stochastic analysis
- `6eq_ssa.c`: simulate the Markov chain model of bipartite virus model using Gillespie's algorithm
- `6eq_rk4_bipartite_var_gamma_3d_tozero.c`: find the initial condition that goes to thge origin projected three dimensional phase plane
- `6eq_rk4_bipartite_var_gamma_R1R2_tozero.c`: find the initial condition that goes to thge origin projected two dimensional phase plane
- `6eq_plot_gillespie.R`: plot the simulation results from `6eq_ssa.c`
- `6eq_ssa_prob_extinction_plot.R`: plot the probability extinction results
- `6eq_plot_errorbars.R`: error bar plot for the probability extinction simulation results

## More info
For more information, do not hesitate to contact me `ruth.kristianingsih30@gmail.com`.


