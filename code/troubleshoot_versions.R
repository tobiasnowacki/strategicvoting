# Load data
library(tidyverse)
library(tictoc)
library(here)
source(here("code/full_header.R"))  # fn's and data
source(here("code/prep_cses.R"))    # data prep
source(here("code/utils/ternary_functions.R"))    # ternary w/o 
source("code/utils/new_sv.R")

# Load Swedish results as is
tic()
results = many_iterations_until_convergence(
                        big_list_na_omit[[18]], 
                        big_list_na_omit[[18]]$v_vec, 
                        0.05, 
                        85, 
                        0.001, 
                        50)
toc()
# Plot v_vec prevalence.
d_by_iter = results[[1]] %>% 
  group_by(iter) %>%
  summarise(prev = mean(tau_plur > 0, na.rm = TRUE))

results[[1]] %>% filter(iter == 12, tau_ > 0, sin_rcv == opt_rcv) %>% head

plot(d_by_iter$iter, d_by_iter$prev)


# Get rid of 'sincerity adjustment' and see if results change.