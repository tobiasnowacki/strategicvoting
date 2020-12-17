# Load data
library(tidyverse)
library(here)
source(here("code/full_header.R"))  # fn's and data
source(here("code/prep_cses.R"))    # data prep
source(here("code/utils/ternary_functions.R"))    # ternary w/o 

# Load Swedish results as is
tic()
results = many_iterations_until_convergence(
                        big_list_na_omit[[32]], 
                        big_list_na_omit[[32]]$v_vec, 
                        0.05, 
                        85, 
                        0.001, 
                        100, ae = FALSE)
toc()
# Plot v_vec prevalence.
d_by_iter = results[[1]] %>% 
  group_by(iter) %>%
  summarise(prev = mean(tau_rcv > 0, na.rm = TRUE))

results[[1]] %>% filter(iter == 12, tau_rcv > 0, sin_rcv == opt_rcv) %>% head

plot(d_by_iter$iter, d_by_iter$prev)


# Get rid of 'sincerity adjustment' and see if results change.