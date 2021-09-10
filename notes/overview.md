# Order of scripts to be run

server.sh

* produces iterated simulations
* needs to be run with all parameter values

all_gather.sh

* produces summary plot, v.vec distances plots and v.vec ternary plot

random_iter.sh

* produces iterated simulations with random starting points
* needs to be run with seeds 1 to 10

random_iter_process.sh

* produces quantile plot and selected cases v.vec ternary plot for random starting locations

after_gather.sh

  - new_convergence.R
  * produces IRV convergence plots 
  * needs to be run once after all_gather.sh
  
  - new_summary.R
  * produces summary plot grid for main paper
  * needs to be run once after all_gather.sh

conjecture.sh

* produces conjecture test(s) for baseline case
