### -------------------------
### LOAD DEPENDENCIES -------

# install.packages("here")
library(here)                       # to get dir
source(here("code/full_header.R"))  # fn's and data
source(here("code/prep_cses.R"))    # data prep

# Set values
lambda = 0.05
epsilon_thresh <- 0.001                      # threshold (epsilon)
s_val = 85
max_iter_val = 160

# Run 2nd case to figure out what's wrong
out = many_iterations_until_convergence(
                        big_list_na_omit[[2]], 
                        big_list_na_omit[[2]]$v_vec, 
                        lambda, 
                        s_val, 
                        epsilon_thresh, 
                        160, ae = TRUE)

# ok, so the issue appears at the 60th iteration. why?
out_59 = one_iteration(big_list_na_omit[[2]], as.numeric(out$rcv_v_vec[59, ]), lambda, s_val, ae)

out_60 = one_iteration(big_list_na_omit[[2]], as.numeric(out$rcv_v_vec[60, ]), lambda, s_val, ae)


test_out = sv(big_list_na_omit[[2]]$U, 
   big_list_na_omit[[2]]$weights, 
   out_60$rcv_vec, 
   85,
   "AV")

glimpse(test_out)

irv_election(n = 1) %>%
    election_event_probs(method = "mc",
                         num_sims = 400000,
                         mc_method = "density",
                         alpha = rep(1/6, times = 6))


irv <- irv_election(n = 10000) 
alpha6 <- c(.3, .05, .2, .15, .1, .2)*85
irv %>% 
  election_event_probs(method = "mc", alpha = alpha6) -> pps

           pps = lapply(pps, function(x) {
            x$integral = x$integral * 10000
            return(x)})
pps


ß
out_60 = one_iteration(big_list_na_omit[[2]], as.numeric(out$rcv_v_vec[61, ]), lambda, s_val, ae)

out_60_vec = out_60$rcv_vec
out_60_br = out_60$rcv_best_response

glimpse(out_60)
out_60_vec


one_iteration(big_list_na_omit[[2]], )