library(tidyverse)
library(pivotprobs)
library(gtools)

# Load data
load("output/big_list_4_parties.RData")
source("code/prep_cses.R") # data prep
big_list_na_omit_4 <- big_list_na_omit
load("output/big_list_2.RData")
source("code/prep_cses.R") # data prep
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")

source("code/utils/new_sv.R")
source("code/utils/sv_four.R")

big_list_na_omit[[5]]$v_vec %>% sum

v.vec <- rdirichlet(1, c(1, 1, 1, 1, 1, 1))
s <- 85

pps <- irv_election(n = 1) %>%
    election_event_probs(method = "en", alpha = (v.vec * s)) 

pmat <- pps %>% combine_P_matrices

v.vec.4 <- rdirichlet(1, rep(1, 24))
mc_sims <- simulate_ordinal_results_from_dirichlet(
    k = 4,
    n = 20000,
    alpha = (v.vec.4 * s)
)

pps4 <- mc_sims %>% irv_pivot_probs_four_cands(n = 1, reporting = 0)
pmat4 <- pps4 %>% combine_P_matrices()

pmat4
pmat

out <- sv(
    U = big_list_na_omit[[5]]$U, 
    s = 85,
    rule = "AV"
)

out_four <- sv_four(
    U = big_list_na_omit_4[[5]]$U, 
    s = 85,
    rule = "AV"
)
