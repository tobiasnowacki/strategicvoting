# Simulation to assess bias from 'fixes'

library(tidyverse)
library(pivotprobs)
library(gtools)

# Simulate v.vecs
vvecdf <- rdirichlet(n = 100000, c(1, 1, 1, 1) )
s <- 85

head(vvecdf)

vvec_test <- c(0.627, 0.365, 0.00467, 0.002287)

pps <- plurality_election(k = 4, n = 1000) %>%
    election_event_probs(
        method = "sc",
        alpha = (vvec_test * s),
        drop_dimension = TRUE,
        merge_adjacent_pivot_events = TRUE
    ) 

vvec_test_amended <- vvec_test + rep(0.002, 4)

pps <- plurality_election(k = 4, n = 1000) %>%
    election_event_probs(
        method = "sc",
        alpha = (vvec_test_amended * s),
        drop_dimension = TRUE,
        merge_adjacent_pivot_events = TRUE
    ) 

int1 <- pps %>% map_dbl("integral")

pps_ev <- plurality_election(k = 4, n = 1000) %>%
    election_event_probs(
        method = "ev",
        alpha = (vvec_test * s)
    ) 

int2 <- pps_ev %>% map_dbl("integral")

pps_sim <- plurality_election(k = 4, n = 1000) %>%
    election_event_probs(
        method = "mc",
        num_sims = 10500000,
        alpha = (vvec_test * s)
    ) 

int3 <- pps_sim %>% map_dbl("integral")

rbind(int1, int2, int3)

glimpse(pps)
