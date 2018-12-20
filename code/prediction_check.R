##################################################
## Project: Strategic Voting in RCV
## Script purpose: Checking theoretical predictions
## Date: 19/12/2018
## Author:
##################################################

### 
### Dependencies
###

library(here)
library(gtools)
library(ggtern)
library(stargazer)

source(here("utils/functions.r"))
source(here("utils", "av_pivotal_probs_analytical_general_v2.r"))
sim_appr2 <- here("utils", "general_iteration_simulation_approach.r")
source(sim_appr2)
sv_file <-  here("utils/sv.r")
source(sv_file)
source(here("utils", "plot_results_on_ternary_mirror.r"))

load(here("..", "output", "cses_big_list_2.RData"))

### 
### Analysis
###

# Create list with v_vecs from CSES utility dfs:
big_list_na_omit <- lapply(big_list, function(x) remove_nas(x))
sin_vote_list <- lapply(big_list_na_omit, function(x) sincere.vote.mat.from.U(x$U, rule = "AV"))
v_vec_list <- list()
for(i in 1:length(sin_vote_list)){
  weights <- big_list_na_omit[[i]]$weights
  v_vec <- ballot.props.from.vote.mat.and.weights(sin_vote_list[[i]], weights)
  big_list_na_omit[[i]]$v_vec <- as.numeric(v_vec)
}

# Drop NA cases
names(big_list)[[39]]
names(big_list)[[145]]
big_list_na_omit[[39]] <- NULL
big_list_na_omit[[144]] <- NULL

# Set level of uncertainty
s <- 85

# For each case:
# 1. get v.vec 
# 1. get "transition matrix"
# 2. get pivotal probabilities

results <- list()
for(i in 1:length(big_list_na_omit)){
	print(i)
	this_list <- big_list_na_omit[[i]]
	pprobs <- av.pivotal.event.probs.general(
		c(this_list$v_vec, 0, 0, 0),
		rep(s, 4)
		)
	class <- classify.vec(this_list$v_vec)
	sv_obj <- convert_andy_to_sv_item_two(
		this_list$U, 
		this_list$weights, 
		s, 
		this_list$v_vec)
	trans_matrix <- vote_matrix_weighted(
		sv_obj,
		type = "rcv", 
		weights = big_list_na_omit[[i]]$weights
		)
	results[[i]] <- list(this_list$v_vec, pprobs, sv_obj, trans_matrix, class)
}

# Classification

# Create matrix of v_vecs
v_vec_df <- t(sapply(results, `[[`, 1))

get_strat_vote_vec <- function(vote_mat){
	vote_mat <- t(apply(vote_mat, 1, function(x) x / sum(x)))
	return(c(vote_mat[1, 3], vote_mat[1, 5], vote_mat[2, 3], vote_mat[2, 5], 
		vote_mat[3, 1], vote_mat[3, 6], vote_mat[4, 1], vote_mat[4, 6], 
		vote_mat[5, 2], vote_mat[5, 4], vote_mat[6, 2], vote_mat[6, 4]))
}

v_vec_df <- as.data.frame(cbind(v_vec_df, t(sapply(results, function(x) get_strat_vote_vec(x[[4]])))))
names(v_vec_df) <- c("ABC", "ACB", "BAC", "BCA", "CAB", "CBA", 
	"ABCBAC", "ABCCAB", "ACBBAC", "ACBCAB", 
	"BACABC", "BACCBA", "BCAABC", "BCACBA",
	"CABACB", "CABBCA", "CBAACB", "CBABCA")

v_vec_df$type <- sapply(results, `[[`, 5)


ggplot(v_vec_df[v_vec_df$type == "SP(B)", ]) +
	geom_density(aes(CBABCA)) + 
	geom_density(aes(ABCCAB), lty = "dotted") + 
	geom_density(aes(ACBCAB), lty = "dashed")

ggplot(v_vec_df[v_vec_df$type == "SP(C)", ]) +
	geom_density(aes(BCACBA)) + 
	geom_density(aes(ABCBAC), lty = "dotted") + 
	geom_density(aes(ACBBAC), lty = "dashed")	

ggplot(v_vec_df[v_vec_df$type == "DM(B)", ]) +
	geom_density(aes(BCACBA)) + 
	geom_density(aes(ABCBAC), lty = "dotted") + 
	geom_density(aes(ACBBAC), lty = "dashed")

ggplot(v_vec_df[v_vec_df$type == "SP(A)", ]) +
	geom_density(aes(BACCBA)) + 
	geom_density(aes(BCACBA), lty = "dotted") + 
	geom_density(aes(CABACB), lty = "dashed")

ggplot(v_vec_df[v_vec_df$type == "DM(A)", ]) +
	geom_density(aes(BACCBA)) + 
	geom_density(aes(BCACBA), lty = "dotted") + 
	geom_density(aes(CABACB), lty = "dashed")

hist(v_vec_df[v_vec_df$type == "SP(B)", ]$BAC / (v_vec_df[v_vec_df$type == "SP(B)", ]$BAC + v_vec_df[v_vec_df$type == "SP(B)", ]$BCA))