library(tidyverse)
source("code/archive/archive_utils/av_pivotal_probs_analytical_general_v2.r")
source("code/archive/archive_utils/plurality_pivotal_probabilities_analytical.r")


setwd("~/Dropbox/research/strategic_voting/code")
source("~/Dropbox/research/strategic_voting/code/utils/EU_given_piv_probs_and_utility.r")
source("~/Dropbox/research/strategic_voting/code/utils/av_pivotal_probs_analytical_general_v2.r")
source("~/Dropbox/research/strategic_voting/code/utils/plurality_pivotal_probabilities_analytical.r")
source("~/Dropbox/research/strategic_voting/code/utils/general_iteration_simulation_approach.r")


this.U = big_list[[32]]$U %>% as.matrix
this.weight = big_list[[32]]$weights
the.floor = 0
out = iterated_best_response_sequence(U = big_list_na_omit[[32]]$U %>% as.matrix, 
                                s = 85, 
                                weights = big_list_na_omit[[32]]$weights, 
                                rule = "AV", 
                                lambda = 0.05, 
                                until.convergence = F, 
                                max.iterations = 150, sincere.proportion = 0, candidates = c("a", "b", "c"), ballots = c("abc", "acb", "bac", "bca", "cab", "cba"), the.floor = 0, noisy = F)

test.vvec = c(results[[2]][12, ] %>% unlist)
test = sv(big_list_na_omit[[32]]$U, big_list_na_omit[[32]]$weights, v.vec = test.vvec, 85, "AV") 
test2 = sv(big_list_na_omit[[32]]$U, big_list_na_omit[[32]]$weights, v.vec = test2.vvec, 85, "AV") 
test$piv.probs
test2$piv.probs
P_mat_at_pivotal_events(test$piv.probs, rule = "AV")
test$eu.mat[585, ]
colSums(test$V.mat)

[1] "best.response.v.vec" "sincere.vote.mat"    "eu.by.ballot"
[4] "weights"             "V.mat"               "P.mat"
[7] "v.vec.before"        "probability.pivotal" "distance.from.last"
out[[13]][[2]][776, ]
out[[13]][[3]][776, ]
out[[13]][[5]][776, ]
out[[13]][[6]]

out[[11]][[3]][585, ]
out[[11]][[5]][585, ]
out[[12]][[5]] %>% colSums
out[[11]][[6]]
out[[11]][[1]]

# test2.vvec = 0.95 * test.vvec + 0.05 * out[[11]][[1]]
test2.vvec = out[[11]][7] %>% unlist
out[[23]][7] %>% unlist
c(results[[2]][23, ] %>% unlist)

mean(test$tau > 0)
mean(test$opt.votes.strategic != test$opt.votes.sincere)

expected.benefit.mat.from.ibrs(out[[11]])

out_eb = expected.benefit.mat.from.ibrs(out)
out_eb[, 13]
out_eb[, 11][585]
test$tau[780] > 0
test$tau[780]
mean(test$tau > 0, na.rm = TRUE)
test$tau
apply(out_eb, 2, function(x) mean(x > 0, na.rm = TRUE))
glimpse(out_eb)